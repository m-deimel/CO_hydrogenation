import re
from copy import deepcopy
from string import ascii_letters

import numpy as np


class BEPProcessHolder(object):
    """
    Serves as a container for processes. Manipulate these processes during the
    building of the model to include nearest neighbors as bystanders and
    modify the otf_rate accordingly.
    """
    error_messages = {
        'flag_self_consistency': 'Self consistency error for bystanders in process: %s',
    }
    warning_messages = {
        'missing_alpha': 'Warning: Missing alpha for %s. Setting alpha to 0.5',
    }

    def __init__(self):
        """ save basic processes and bystanders of sites """
        self.processes = []
        self.site_bystanders = {}
        self.param_list = []
        self.interaction_energy_pattern = (
            r'I_([A-Z0-9]+)_{bystander_site}_{species}_{species_site}'
            '|I_{species}_{species_site}_([A-Z0-9]+)_{bystander_site}'
        )
        self.self_interaction_energy_pattern = r'I_{c1_species}_{c1_site}_{c2_species}_{c2_site}'

    def add_process(self, **kwargs):
        """ add basic process to this holder """
        self.processes.append(kwargs)

    def add_param_list(self, param_list):
        self.parameter_list = param_list

    def add_site_bystanders(self, site_name, site, bystanders):
        """ adds a list of bystanders for a given site """
        self.site_bystanders[site_name] = [site, bystanders]

    def get_site_bystanders(self, coord, flag):
        """
        generate bystanders of the given coordinate,
        by using the name and offset of self.site_bystanders
        """
        bystanders = deepcopy(self.site_bystanders[coord.name][1])
        offset = coord.offset - self.site_bystanders[coord.name][0].offset
        for bystander in bystanders:
            bystander.coord.offset += offset
            bystander.flag += '_' + flag
        return bystanders

    def get_rate_pairs(self, sign, bystander_list, species, coord):
        """
        Given a bystander list and a species at a given coordinate, we
        calculate the rate modification summands.
        A rate modification summand is generated, if a species at a site
        has an interaction to a bystander in bystander list.
        We use the kmos project parameter list pt.parameter_list, to lookup
        these interactions.

        Parameters:
        -----------
        sign: '+' or '-'
            '+' for initial state and '-' for final state
        bystander_list: 1-D list of kmos bystanders
        species: kmos species
        coord: kmos coord

        Returns:
        --------
        rate_pairs: List of strings to build a summand of the otf rate expression
        bystander_list: List of bystanders with adjusted allowed species
        """
        rate_pairs = []
        for bystander in bystander_list:
            # species and initial_coordinate given
            matched_species = []
            pattern = self.interaction_energy_pattern.format(
                bystander_site=bystander.coord.name,
                species=species,
                species_site=coord.name
            )

            for parameter in self.parameter_list:
                matches = re.match(pattern, parameter.name)
                if matches:
                    group1, group2 = matches.groups()
                    if group1:
                        matched_species.append(group1)
                    else:
                        matched_species.append(group2)

                    # summands within reaction rate adjustment, e.g. lines 2 and 3 in the following comment
                    # otf_rate=('base_rate*exp(beta*alpha_CO_des*('
                    #           'I_CO_cus_CO_cus*nr_CO_nn_cus + I_O_cus_CO_cus*nr_O_nn_cus +'
                    #           'I_CO_cus_CO_br*nr_CO_nn_br + I_O_br_CO_cus*nr_O_nn_br'
                    #           ')*eV)')
                    # rate_pairs.append(sign + parameter.name + '*nr_' + matched_species[-1] + '_' + bystander.flag)
                    rate_pairs.append((sign, parameter.name, matched_species[-1], bystander.flag))
            if bystander.allowed_species:
                matched_species = list(set(bystander.allowed_species + matched_species))
            bystander.allowed_species = matched_species
        return rate_pairs, bystander_list

    @staticmethod
    def rate_modification(rate_pairs):
        """
        construct otf_rate modification by interaction/species/site pairs

        Parameters:
        -----------
        rate_pairs: 1-D list of 4-tuples
            [(sign, Interaction, bystander species, bystander flag), ]

        Returns:
        --------
        rate_modification: str
            string concatenation of all bystander summands. One summand looks like:
            "[sign][Interaction]*nr_[bystander species]_[bystander flag]"
        """
        # remove full duplicates
        rate_pairs = list(set(rate_pairs))

        rate_modification = ''
        for sign, interaction, species, flag in rate_pairs:
            rate_modification += sign + interaction + '*nr_' + species + '_' + flag
        return rate_modification

    @staticmethod
    def reduce_rate_modification(rate_modification):
        """
        Reduce the rate_modification arithmetically as much as possible.

        Usually the rate_modification looks something like:
        "Interaction1_Energy*Count_of_Interactions1 + Interaction2_Energy*Count_of_Interactions2 + ..."
        Often there are summands which cancel out completely (usually the case for diffusions).
        Sometimes we can apply a distributive law like:
        "Interaction1_Energy*Count_of_Interactions1_initial_state -
         Interaction1_Energy*Count_of_Interactions1_final_state =
         Interaction1_Energy*(Count_of_Interaction1_initial_state-Count_of_Interaction1_final_state)"
        Sometimes there are "self-interactions" between reactants / products which are represented by
        "Interaction1_Energy" - we can apply the distributive law here, too, via
        "Interaction1_Energy*1"

        Parameters:
        -----------
        rate_modification: str

        Returns:
        --------
        reduced_rate_modification: str
        """
        rate_modification = rate_modification.strip().replace(" ", "")

        # explicitly insert a '+' sign for first term for more robust term matching
        if not (rate_modification.startswith('+') or rate_modification.startswith('-')):
            rate_modification = '+' + rate_modification

        # ['+a*b', '-a*c', '+a'] = re.findall(term_pattern, '+a*b-a*c+a')
        term_pattern = r'[\+|-]?[a-zA-Z0-9_\*]+'
        terms = re.findall(term_pattern, rate_modification)

        # for each factor in each, create a list of corresponding factors
        # a -> [+b, -c, +1]
        # b -> [+a]
        # c -> [-a]
        # Use this as heuristic for the shortest possible rate_modification
        gathered_factors = {}
        for term in terms:
            sign = term[0]
            if '*' in term[1:]:
                factor1, factor2 = term[1:].split('*')
            else:
                factor1, factor2 = term[1:], '1'

            if factor1 not in gathered_factors:
                gathered_factors[factor1] = []
            gathered_factors[factor1].append(sign+factor2)

            if factor2 != '1':
                if factor2 not in gathered_factors:
                    gathered_factors[factor2] = []
                gathered_factors[factor2].append(sign+factor1)

        # The following steps could be done more computationally efficient
        # using trees, however computational efficiency is no problem at all here

        # remove factors, which cancel each other
        for key, factors in gathered_factors.items():
            to_remove = []
            for i, factor1 in enumerate(factors):
                sign1 = factor1[0]
                tmp_factor1 = factor1[1:]
                for factor2 in factors[i+1:]:
                    sign2 = factor2[0]
                    tmp_factor2 = factor2[1:]
                    if tmp_factor1 == tmp_factor2 and sign1 != sign2:
                        to_remove.extend([factor1, factor2])
            for factor in to_remove:
                factors.remove(factor)
            if not factors:
                gathered_factors.pop(key)

        # now always take the key with biggest length and build
        # the reduced rate expression
        reduced_rate_modification = ''
        while gathered_factors:
            max_length = -1
            max_key = None
            for key, value in gathered_factors.items():
                if len(value) > max_length:
                    max_length = len(value)
                    max_key = key

            key = max_key
            factors = gathered_factors.pop(key)

            if max_length == 1:
                reduced_rate_modification += factors[0][0] + key + '*' + factors[0][1:]
            else:
                reduced_rate_modification += '+' + key + '*(' + ''.join(factors) + ')'

            # update gathered factors
            for factor in factors:
                sign = factor[0]
                tmp_factor = factor[1:]
                if tmp_factor != '1':
                    gathered_factors[tmp_factor].remove(sign+key)
                    if not gathered_factors[tmp_factor]:
                        gathered_factors.pop(tmp_factor)

        return reduced_rate_modification

    def self_interactions(self, sign, condition_list):
        """
        finds interactions between reactants / products and returns
        appropriate interaction energies.

        Parameters:
        -----------
        condition_list: 1-D list of kmos condition

        Returns:
        --------
        interaction_parameters: [I_condition1_condition2, I_condition1_condition3, ..., I_condition2_condition3, ...]
        """
        interaction_parameters = []

        for i, condition1 in enumerate(condition_list):
            for j in range(i+1, len(condition_list)):
                condition2 = condition_list[j]
                parameter_name1 = self.self_interaction_energy_pattern.format(
                    c1_species=condition1.species,
                    c1_site=condition1.coord.name,
                    c2_species=condition2.species,
                    c2_site=condition2.coord.name,
                )
                parameter_name2 = self.self_interaction_energy_pattern.format(
                    c2_species=condition1.species,
                    c2_site=condition1.coord.name,
                    c1_species=condition2.species,
                    c1_site=condition2.coord.name,
                )

                for param in self.parameter_list:
                    if param.name == parameter_name1 or param.name == parameter_name2:
                        interaction_parameters.append(sign + param.name)

        return interaction_parameters

    @staticmethod
    def flatten_bystander_list(bystander_list, process_name):
        """
        Takes a 2d list of bystanders for each involved coordinate.
        If the coordinate name and offset are the same, we can combine the
        allowed species. If a bystander doesn't have any allowed species,
        we may simply drop it.

        Parameters:
        -----------
        bystander_list: [[bystanders_of_coord1], [bystanders_of_coord2], ...]

        Returns:
        --------
        bystander_list: 1-D list of kmos bystanders
        """
        tmp = []
        for bystanders in bystander_list:
            for bystander in bystanders:
                # if there are no allowed species, we can just drop the bystander
                if not bystander.allowed_species:
                    continue

                tmp_index = -1
                for i, item in enumerate(tmp):
                    # if the bystander is already in tmp,
                    # add allowed species to existing bystander
                    if item.coord.name == bystander.coord.name and np.all(item.coord.offset == bystander.coord.offset):
                        tmp_index = i
                        item.allowed_species = list(set(item.allowed_species + bystander.allowed_species))

                        # self consistency check - if it's the same bystander
                        # both item and bystander should have the same flag
                        assert item.flag == bystander.flag, \
                            BEPProcessHolder.error_messages['flag_self_consistency'] % process_name

                if tmp_index == -1:
                    tmp.append(bystander)
        return tmp

    def get_coord_list_bystanders(self, flag_base, coordinates, species):
        """
        Generates a list of bystanders for each coordinate in coordinates.

        Parameters:
        -----------
        flag_base: str
            'initial' or 'final' depending on initial or final state of coordinates
        coordinates: 1-D list of kmos coordinates
        species: 1-D list of kmos species
            [species@coord1, species@coord2, ...]

        Returns:
        --------
        bystanders: 2-D List
            [[bystanders@coord1], [bystanders@coord2], ...]
        """
        bystanders = []
        for i, (coord, species) in enumerate(zip(coordinates, species)):
            if species == 'empty':
                bystanders.append([])
            else:
                bystanders.append(self.get_site_bystanders(coord, flag_base + str(i)))
        return bystanders

    @staticmethod
    def hashable_coordinate(coord):
        """
        To represent coordinates in hashmaps, we require a hashable
        representation of the coordinate. Generate a string out of name
        and offset of the coordinate.

        Parameters:
        -----------
        coord: kmos coordinate

        Returns:
        --------
        hashable_representation: str
        """
        return coord.name + ';'.join(map(str, coord.offset))

    @staticmethod
    def bystander_in_coordinate_list(bystander, coords):
        """
        Helper for filter function.

        Parameters:
        -----------
        bystander: kmos bystander
        coords: 1-D list of kmos coordinates

        Returns:
        --------
        in_list: bool
        """
        in_list = False
        for coord in coords:
            if coord.name == bystander.coord.name and np.all(coord.offset == bystander.coord.offset):
                in_list = True
                break
        return in_list

    def react_otf(self, process, alpha):
        """
        First find all bystanders for each site involved in the given process.
        Given the species and calculated bystanders, the otf rate expression is
        calculated. The alpha for the BEP-relation has to be in the project parameters.

        Paramaters:
        -----------
        process: kmos process
        alpha: str
            Name of the param describing the BEP slope for the given process

        Returns:
        --------
        bystanders: 1-D list of bystanders
        rate_modification: str otf rate
        """

        initial_coordinates = [condition.coord for condition in process['condition_list']]
        initial_species = [condition.species for condition in process['condition_list']]
        final_coordinates = [action.coord for action in process['action_list']]
        final_species = [action.species for action in process['action_list']]

        initial_bystanders = self.get_coord_list_bystanders('initial', initial_coordinates, initial_species)
        final_bystanders = self.get_coord_list_bystanders('final', final_coordinates, final_species)

        # exclude condition / action coordinates from bystanders
        # as we know the species at these positions already
        # (also it is a requirement imposed by the otf backend)
        for x_bystanders in [initial_bystanders, final_bystanders]:
            for i, bystanders in enumerate(x_bystanders):
                x_bystanders[i] = [
                    bystander for bystander in bystanders
                    if not self.__class__.bystander_in_coordinate_list(bystander, initial_coordinates)
                ]

        # find bystanders which are bystanders of multiple conditions / actions
        overlap = {}
        for flag, x_bystanders in [('initial', initial_bystanders), ('final', final_bystanders)]:
            for i, bystanders in enumerate(x_bystanders):
                for bystander in bystanders:
                    coord = self.__class__.hashable_coordinate(bystander.coord)
                    if coord not in overlap:
                        overlap[coord] = {'initial': [], 'final': []}
                    overlap[coord][flag].append(str(i))

        # generate a flag for each bystander to represent
        # all conditions / actions in the flag name
        for key, values in overlap.items():
            flag = key[0]
            if values['initial']:
                flag += '_initial' + ''.join(values['initial'])
            if values['final']:
                flag += '_final' + ''.join(values['final'])
            overlap[key]['flag'] = flag

        # adjust the bystander flags
        # create flags which are not part of other flags...
        # otherwise otf will not generate proper fortran files
        # should be fixed in kmos at some point...
        base = ascii_letters
        unique_flags = {}
        for x_bystanders in [initial_bystanders, final_bystanders]:
            for bystanders in x_bystanders:
                for bystander in bystanders:
                    flag = overlap[self.__class__.hashable_coordinate(bystander.coord)]['flag']
                    if flag not in unique_flags:
                        unique_flags[flag] = base[len(unique_flags)]
                    bystander.flag = unique_flags[flag]

        # calculate rate pairs and adjust allowed_species for each bystander
        # based on the existing lateral interaction energies in the project parameters
        rate_pairs = []
        rate_pair_parameters = [
            ('+', initial_bystanders, initial_species, initial_coordinates),
            ('-', final_bystanders, final_species, final_coordinates),
        ]
        for sign, x_bystanders, x_species, x_coordinates in rate_pair_parameters:
            for bystanders, species, coord in zip(x_bystanders, x_species, x_coordinates):
                tmp, bystanders = self.get_rate_pairs(sign, bystanders, species, coord)
                rate_pairs.extend(tmp)

        # create a flat bystander list, required for the otf framework
        bystanders = self.__class__.flatten_bystander_list(initial_bystanders + final_bystanders, process['name'])
        self_interaction_parameters = self.self_interactions('+', process['condition_list'])
        self_interaction_parameters += self.self_interactions('-', process['action_list'])

        rate_modification = self.__class__.rate_modification(rate_pairs) + ''.join(self_interaction_parameters)
        rate_modification = self.__class__.reduce_rate_modification(rate_modification)

        # after reducing the rate_modification it is possible, that we can drop bystanders
        # or allowed species within these bystanders
        # (often the case for diffusion processes -> canceling elements)
        bystanders = [bystander for bystander in bystanders if '_' + bystander.flag in rate_modification]

        return bystanders, 'base_rate*exp(%s*beta*(%s)*eV)' % (alpha, rate_modification)

    def get_numeric_alpha_value(self, alpha_name, pt):
        """
        Searches self.param_list for alpha and returns it's numeric value.
        """
        for param in self.parameter_list:
            if param.name == alpha_name:
                return float(param.value)

    def add_project_processes(self, pt):
        """ add all processes to the project including otf_rate and bystander_list"""
        for process in self.processes:
            alpha = process.pop('alpha', None)
            if not alpha:
                if 'diff' in process['name']:
                    # BEP slope for diffusion processes of the same site is 0.5
                    # if process['condition_list'][0].coord.name == process['condition_list'][1].coord.name:
                    #     alpha = 'alpha'

                    # for now our assumption is, that diffusion processes always have a BEP slope of alpha=0.5
                    alpha = 'alpha'
                if not alpha:
                    output = ''
                    for condition in process['condition_list']:
                        output += condition.species + '@' + condition.coord.name + ' + '
                    output = output[:-2]
                    output += ' -> '
                    for condition in process['action_list']:
                        output += condition.species + '@' + condition.coord.name + ' + '
                    output = output[:-2]
                    print(BEPProcessHolder.warning_messages['missing_alpha'] % process['name'] + ' - ' + output)

                    alpha = 'alpha'

            numeric_alpha = None
            if alpha.startswith('rev_'):
                numeric_alpha = 1 - self.get_numeric_alpha_value(alpha[4:], pt)
                alpha = '(1-%s)' % alpha[4:]
            else:
                numeric_alpha = self.get_numeric_alpha_value(alpha, pt)

            bystander_list, otf_rate = None, None
            if process['rate_constant'] != '0.0' and abs(numeric_alpha) > 1e-10:
                bystander_list, otf_rate = self.react_otf(process, alpha)

            if bystander_list and otf_rate:
                process['bystander_list'] = bystander_list
                process['otf_rate'] = otf_rate

            pt.add_process(**process)
        return pt


