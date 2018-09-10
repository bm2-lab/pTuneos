#!/usr/bin/env python

'''
Created on 09.21.2015
'''

from __future__ import print_function

from collections import namedtuple
import os
import sqlite3 as lite

DB_FILENAME = 'allele.sqlite'
DB_FULL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), DB_FILENAME)

def open_and_close_sqlite3_connection(func):
    """ | *brief*: Decorator which wraps *func* so if a 'connection' kwarg isn't provided
        |    it creates and passes one.  If one is created, it is closed.
        | *author*: Jivan
        | *created*: 2015-11-24
    """
    def wrapped_func(*args, **kwargs):
        connection = kwargs.get('connection')
        if connection:
            close_connection = False
        else:
            connection = lite.connect(DB_FULL_PATH)
            close_connection = True

        try:
            kwargs.update({'connection': connection})
            ret = func(*args, **kwargs)
        finally:
            if close_connection:
                connection.close()
        return ret
    return wrapped_func


class MHCIAlleleData(object):
    @open_and_close_sqlite3_connection
    def get_method_names(self, allele_name=None, binding_length=None,
                         include_hidden=False, connection=None):
        """ @brief: Returns a list of valid method names for MHCI.
            @author: Jivan
            @since: 2015-10-06
            Some methods are hidden from users, but can still be found by setting
            include_hiddent to True.  As of 2015-12-15, this was just method 'arb'.
        """
        sql = '''
            SELECT DISTINCT mm.name
            FROM mhci_allele ma JOIN mhci_method mm ON ma.method_id = mm.id
            {}
            ORDER BY mm.name;
        '''
        # --- Build a parameterized WHERE clause based on the optional keyword arguments.
        # SQLite uses integers for booleans.  If include_hidden is not false, don't limit by
        #    this parameter.
        include_hidden = 0 if not include_hidden else None
        where_pieces = ['ma.name = ?', 'ma.length = ?', 'mm.hidden = ?']
        where_args = [allele_name, binding_length, include_hidden]
        assert(len(where_pieces) == len(where_args))
        # Remove the corresonding where pieces & where args if a where arg is None.
        filtered_where_pieces = [ wp for wp, wa in zip(where_pieces, where_args) if wa is not None ]
        filtered_where_args = [ wa for wa in where_args if wa is not None ]
        where_clause_content = ' AND '.join(filtered_where_pieces)

        if filtered_where_args:
            where_clause = 'WHERE {}'.format(where_clause_content)
        else:
            where_clause = ''

        # --- Execute the query and return the results.
        sql = sql.format(where_clause)
        cur = connection.cursor()
        cur.execute(sql, filtered_where_args)
        results = cur.fetchall()
        method_names = [row[0] for row in results]
        return method_names

    @open_and_close_sqlite3_connection
    def get_all_allele_names(self, connection=None):
        """ @brief: Returns an ordered list of all allele names.
            @author: Jivan
            @since: 2015-11-05
        """
        sql = 'SELECT DISTINCT name FROM mhci_allele ORDER BY name;'
        cur = connection.cursor()
        cur.execute(sql)
        results = cur.fetchall()
        allele_names = [ row[0] for row in results ]
        return allele_names

    @open_and_close_sqlite3_connection
    def get_allele_names_for_method(self, method_name, binding_length=None, connection=None):
        """ @brief: Returns a list of allele names valid for use with the method indicated
                by \a method_name.
            This is allows get_allowed_peptide_lengths() to validate method_name /
            allele_name combinations.
            @author: Jivan
            @since: 2015-10-08
        """
        # Make method_name's passed case-insensitive.
        method_name = method_name.lower()

        sql = '''
            SELECT DISTINCT ma.name
            FROM mhci_allele ma JOIN mhci_method mm ON ma.method_id = mm.id
            WHERE {};
        '''
        where_pieces = ['mm.name = ?']
        where_args = [method_name]
        if binding_length is not None:
            where_pieces.append('ma.length = ?')
            where_args.append(binding_length)

        where_clause = ' AND '.join(where_pieces)
        sql = sql.format(where_clause)

        cur = connection.cursor()
        cur.execute(sql, where_args)
        results = cur.fetchall()
        allele_names = [ row[0] for row in results ]
        return allele_names

    @open_and_close_sqlite3_connection
    def get_allele_names(self, species=None, method_name=None, frequency_cutoff=None,
                         connection=None, include_hidden=False):
        """ @brief Returns a list of allele names given a \a species and optionally
                filtered by \a method_name or \a frequency_cutoff.
            @author: Jivan
            @since: 2015-10-06
        """
        # Make method_name's passed case-insensitive.
        if method_name: method_name = method_name.lower()

        if not species:
            raise ValueError('species is a required keyword argument')
        if not species or species not in self.get_species_list():
            raise ValueError('Invalid species: {}'.format(species))
        if method_name and method_name not in self.get_method_names(include_hidden=include_hidden):
            raise ValueError('Invalid method_name: {}'.format(method_name))

        species_clause = 'ma.species = ?'
        method_clause = 'mm.name = ?'
        frequency_clause = 'ma.frequency >= ?'
        sql_where_parts = [species_clause]
        sql_params = [species]
        if method_name:
            sql_where_parts.append(method_clause)
            sql_params.append(method_name)
        if frequency_cutoff:
            sql_where_parts.append(frequency_clause)
            sql_params.append(frequency_cutoff)

        sql_where_clause = ' AND '.join(sql_where_parts)

        sql = '''
            SELECT DISTINCT ma.name
            FROM mhci_allele ma JOIN mhci_method mm ON ma.method_id = mm.id
            WHERE {}
            ORDER BY ma.name;
        '''.format(sql_where_clause)

        cur = connection.cursor()
        cur.execute(sql, sql_params)
        results = cur.fetchall()
        allele_names = [row[0] for row in results]
        return allele_names

    @open_and_close_sqlite3_connection
    def get_allele_frequencies(self, allele_name_list, connection=None):
        """ | *brief*: Returns a list of frequencies for the *allele_name_list* passed.
            | *author*: Jivan
            | *created*: 2015-11-25

            If the allele has no frequency or is invalid it's returned value will be None.
        """
        sql = '''
            SELECT DISTINCT name, frequency
            FROM mhci_allele
            WHERE frequency NOT NULL AND name in ({});
        '''.format(','.join(['?' for i in range(len(allele_name_list))]))
        cur = connection.cursor()
        cur.execute(sql, allele_name_list)
        results = cur.fetchall()
        f_lookup = {row[0]: row[1] for row in results}
        frequencies = [ f_lookup[allele_name] if allele_name in f_lookup else None
                            for allele_name in allele_name_list ]
        return frequencies

    @open_and_close_sqlite3_connection
    def get_species_list(self, method_name=None, connection=None, include_hidden=False):
        """ @brief Returns a list of species for MHCI optionally filtered by \a method_name.
            @author: Jivan
            @since: 2015-10-06
        """
        # Make method_name's passed case-insensitive.
        if method_name: method_name = method_name.lower()

        if method_name and method_name not in self.get_method_names(include_hidden=include_hidden):
            raise Exception('Invalid method name: {}'.format(method_name))

        sql_params = []
        if method_name:
            sql_where_clause = 'WHERE mm.name = ?'
            sql_params.append(method_name)
        else:
            sql_where_clause = ''

        sql = '''
            SELECT DISTINCT ma.species
            FROM mhci_allele ma JOIN mhci_method mm ON ma.method_id = mm.id
            {}
            ORDER BY ma.species
        '''.format(sql_where_clause)

        cur = connection.cursor()
        cur.execute(sql, sql_params)
        results = cur.fetchall()
        species_list = [row[0] for row in results]
        return species_list

    @staticmethod
    def get_species_for_allele_name(allele_name):
        """ @brief: Returns \a species from \a given allele name.
            @author: Dorjee
            @since: 2016-11-01
        """
        # Remove all white-spaces
        allele_name = allele_name.strip()
        # Raise error if allele_name is not provided
        if not allele_name:
            raise ValueError('allele_name is a required keyword argument')

        species_indicator = allele_name[:2].lower()
        species_by_indicator = {
            'hl': 'human',
            'h-': 'mouse',
            'pa': 'chimpanzee',
            'ma': 'macaque',
            'go': 'gorilla',
            'sl': 'pig',
            'bo': 'cow',
            'rt': 'rat',
        }
        species = species_by_indicator[species_indicator] if species_indicator in species_by_indicator else None
        return species

    @open_and_close_sqlite3_connection
    def get_allowed_peptide_lengths(self, method_name, allele_name, connection=None, include_hidden=False):
        """ @brief Returns the valid peptide lengths for the \a species & \a method_name given.
            @author: Jivan
            @since: 2015-10-06
        """
        # Make method_name's passed case-insensitive.
        method_name = method_name.lower()

        if method_name not in self.get_method_names(include_hidden=include_hidden):
            raise ValueError('invalid method_name: {}'.format(method_name))
        valid_allele_names = self.get_allele_names_for_method(method_name)
        if allele_name not in valid_allele_names:
            raise ValueError('invalid allele_name: {}'.format(allele_name))

        sql = '''
            SELECT DISTINCT ma.length
            FROM mhci_allele ma JOIN mhci_method mm ON ma.method_id = mm.id
            WHERE mm.name = ? AND ma.name = ?
            ORDER BY ma.length;
        '''
        cur = connection.cursor()
        cur.execute(sql, [method_name, allele_name])
        results = cur.fetchall()
        lengths = [ row[0] for row in results ]
        return lengths

    def get_reference_set(self):
        """ @brief Returns a list of (<allele name>, <binding length) tuples representing
                the MHCI reference set.
            @author: Jivan
            @since: 2015-10-06
        """
        raw_reference_set = [('HLA-A*01:01', '9'), ('HLA-A*01:01', '10'), ('HLA-A*02:01', '9'), ('HLA-A*02:01', '10'), ('HLA-A*02:03', '9'), ('HLA-A*02:03', '10'), ('HLA-A*02:06', '9'), ('HLA-A*02:06', '10'), ('HLA-A*03:01', '9'), ('HLA-A*03:01', '10'), ('HLA-A*11:01', '9'), ('HLA-A*11:01', '10'), ('HLA-A*23:01', '9'), ('HLA-A*23:01', '10'), ('HLA-A*24:02', '9'), ('HLA-A*24:02', '10'), ('HLA-A*26:01', '9'), ('HLA-A*26:01', '10'), ('HLA-A*30:01', '9'), ('HLA-A*30:01', '10'), ('HLA-A*30:02', '9'), ('HLA-A*30:02', '10'), ('HLA-A*31:01', '9'), ('HLA-A*31:01', '10'), ('HLA-A*32:01', '9'), ('HLA-A*32:01', '10'), ('HLA-A*33:01', '9'), ('HLA-A*33:01', '10'), ('HLA-A*68:01', '9'), ('HLA-A*68:01', '10'), ('HLA-A*68:02', '9'), ('HLA-A*68:02', '10'), ('HLA-B*07:02', '9'), ('HLA-B*07:02', '10'), ('HLA-B*08:01', '9'), ('HLA-B*08:01', '10'), ('HLA-B*15:01', '9'), ('HLA-B*15:01', '10'), ('HLA-B*35:01', '9'), ('HLA-B*35:01', '10'), ('HLA-B*40:01', '9'), ('HLA-B*40:01', '10'), ('HLA-B*44:02', '9'), ('HLA-B*44:02', '10'), ('HLA-B*44:03', '9'), ('HLA-B*44:03', '10'), ('HLA-B*51:01', '9'), ('HLA-B*51:01', '10'), ('HLA-B*53:01', '9'), ('HLA-B*53:01', '10'), ('HLA-B*57:01', '9'), ('HLA-B*57:01', '10'), ('HLA-B*58:01', '9'), ('HLA-B*58:01', '10')]
        AlleleLengthTuple = namedtuple('AlleleLengthTuple', ['allele_name', 'binding_length'])
        # Convert raw reference set to named tuples for easier debugging & use of returned value.
        reference_set = [AlleleLengthTuple(*alt) for alt in raw_reference_set]
        return reference_set


class MHCIIAlleleData(object):

    @open_and_close_sqlite3_connection
    def get_method_names(self, connection=None):
        """ @brief: Returns a list of valid method names for MHC II.
            @author: Sinu
            @since: 2015-10-07
        """
        sql = 'SELECT DISTINCT name FROM mhcii_method ORDER BY name;'
        cur = connection.cursor()
        cur.execute(sql)
        results = cur.fetchall()
        method_names = [row[0] for row in results]
        return method_names

    @open_and_close_sqlite3_connection
    def get_locus_names(self, connection=None):
        """ @brief: Returns a list of valid locus names for MHC II.
            @author: Sinu
            @since: 2015-10-07
        """
        sql = 'SELECT DISTINCT locus FROM mhcii_allele_a ORDER BY locus;'
        cur = connection.cursor()
        cur.execute(sql)
        results = cur.fetchall()
        locus_names = [row[0] for row in results]
        return locus_names

    @open_and_close_sqlite3_connection
    def get_locus_names_for_method(self, method_name=None, connection=None):
        """ @brief: Returns a list of valid locus names given a \a method for MHC II.
            @author: Dorjee
            @since: 2015-03-03
        """
        if not method_name:
            raise ValueError('method_name is a required keyword argument')
        method_name = method_name.lower()

        sql = '''
            SELECT distinct locus 
            FROM mhcii_method mm JOIN mhcii_allele_a maa ON mm.id == maa.method_id 
            WHERE mm.name = ?
            ORDER BY locus;
        '''
        cur = connection.cursor()
        cur.execute(sql, [method_name])
        results = cur.fetchall()
        locus_names = [row[0] for row in results]
        return locus_names

    @open_and_close_sqlite3_connection
    def get_allele_names(self, method_name=None, locus_name=None, connection=None):
        """ @brief: Returns a list of allele names given a \a method and a \a locus.
            @author: Sinu
            @since: 2015-10-07
        """
        if not method_name:
            raise ValueError('method_name is a required keyword argument')
        method_name = method_name.lower()
        if method_name not in self.get_method_names():
            raise ValueError('Invalid method_name: {}'.format(method_name))

        if locus_name and locus_name not in self.get_locus_names():
            raise ValueError('Invalid locus_name: {}'.format(locus_name))

        where_parts = ['mm.name = ?']
        where_args = [method_name]
        if locus_name:
            where_parts.append('ma.locus = ?')
            where_args.append(locus_name)

        where_clause = ' AND '.join(where_parts)
        sql_template = '''
            SELECT DISTINCT ma.name || ifnull("/" || mb.name, '') AS merged_name
            FROM mhcii_method mm
                JOIN mhcii_allele_a ma ON mm.id = ma.method_id
                LEFT JOIN mhcii_allele_b mb ON mb.allele_a_name = ma.name AND mb.method_id = ma.method_id
            WHERE {}
            ORDER BY merged_name;
        '''
        sql = sql_template.format(where_clause)

        cur = connection.cursor()
        cur.execute(sql, where_args)
        results = cur.fetchall()
        allele_names = [row[0] for row in results]
        return allele_names

    @open_and_close_sqlite3_connection
    def get_alpha_chain(self, method_name=None, locus_name=None, connection=None):
        """ @brief: Returns a list of alpha chains given a \a method and a \a locus.
            @author: Dorjee
            @since: 2016-11-03
        """
        if not method_name:
            raise ValueError('method_name is a required keyword argument')
        method_name = method_name.lower()
        if method_name not in self.get_method_names():
            raise ValueError('Invalid method_name: {}'.format(method_name))

        if not locus_name:
            raise ValueError('locus_name is a required keyword argument')
        elif locus_name not in self.get_locus_names():
            raise ValueError('Invalid locus_name: {}'.format(locus_name))

        sql = '''
            SELECT DISTINCT ma.name
            FROM mhcii_allele_a ma JOIN mhcii_method mm ON ma.method_id == mm.id 
            WHERE ma.locus = ? AND mm.name = ?
            ORDER BY ma.name;       
        '''

        cur = connection.cursor()
        cur.execute(sql, [locus_name, method_name])
        results = cur.fetchall()
        allele_names = [row[0] for row in results]
        return allele_names


    @open_and_close_sqlite3_connection
    def get_beta_chain(self, method_name=None, allele_name=None, connection=None):
        """ @brief: Returns a list of beta chains given a \a method and a \a allele.
            @author: Sinu
            @since: 2015-10-07
        """
        if not method_name:
            raise ValueError('method_name is a required keyword argument')
        method_name = method_name.lower()
        if method_name not in self.get_method_names():
            raise ValueError('Invalid method_name: {}'.format(method_name))

        if not allele_name:
            raise ValueError('allele_name is a required keyword argument')

        sql = '''
            SELECT DISTINCT mab.name
            FROM mhcii_allele_b mab JOIN mhcii_method mm ON mab.method_id == mm.id
            WHERE mm.name == ? AND mab.allele_a_name = ?
            ORDER BY mab.name;  
        '''

        cur = connection.cursor()
        cur.execute(sql, [method_name, allele_name])
        results = cur.fetchall()
        allele_names = [row[0] for row in results]
        return allele_names

    @open_and_close_sqlite3_connection
    def get_alpha_chains_for_locus(self, locus_name=None, connection=None):
        """ @brief: Returns a list of beta chains given a \a locus.
            @author: Dorjee
            @since: 2015-02-04
        """
        if not locus_name:
            raise ValueError('locus_name is a required keyword argument')
        locus_name = locus_name.upper()

        sql = 'SELECT DISTINCT name FROM mhcii_allele_a WHERE locus = ?'

        cur = connection.cursor()
        cur.execute(sql, [locus_name])
        results = cur.fetchall()
        allele_names = [row[0] for row in results]
        return allele_names

    @open_and_close_sqlite3_connection
    def get_beta_chains_for_locus(self, locus_name=None, connection=None):
        """ @brief: Returns a list of beta chains given a \a locus.
            @author: Dorjee
            @since: 2015-02-02
        """
        if not locus_name:
            raise ValueError('locus_name is a required keyword argument')
        locus_name = locus_name.upper()

        sql = 'SELECT DISTINCT name FROM mhcii_allele_b WHERE locus = ? ORDER BY name;'

        cur = connection.cursor()
        cur.execute(sql, [locus_name])
        results = cur.fetchall()
        allele_names = [row[0] for row in results]
        return allele_names

    def get_reference_set(self):
        """ @brief Returns a list of allele names representing the MHCII reference set.
            @author: Sinu
            @since: 2015-10-07
        """
        reference_set = ['HLA-DPA1*01/DPB1*04:01', 'HLA-DPA1*01:03/DPB1*02:01', 'HLA-DPA1*02:01/DPB1*01:01', 'HLA-DPA1*02:01/DPB1*05:01', 'HLA-DPA1*03:01/DPB1*04:02', 'HLA-DQA1*01:01/DQB1*05:01', 'HLA-DQA1*01:02/DQB1*06:02', 'HLA-DQA1*03:01/DQB1*03:02', 'HLA-DQA1*04:01/DQB1*04:02', 'HLA-DQA1*05:01/DQB1*02:01', 'HLA-DQA1*05:01/DQB1*03:01', 'HLA-DRB1*01:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*04:05', 'HLA-DRB1*07:01', 'HLA-DRB1*08:02', 'HLA-DRB1*09:01', 'HLA-DRB1*11:01', 'HLA-DRB1*12:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01', 'HLA-DRB3*01:01', 'HLA-DRB3*02:02', 'HLA-DRB4*01:01', 'HLA-DRB5*01:01']
        return reference_set


class MHCNPAlleleData(object):

    @open_and_close_sqlite3_connection
    def get_method_names(self, connection=None):
        """ @brief: Returns a list of valid method names for MHCNP.
            @author: Dorjee
            @since: 2015-10-07
        """
        sql = 'SELECT DISTINCT name FROM mhcnp_method ORDER BY name;'
        cur = connection.cursor()
        cur.execute(sql)
        results = cur.fetchall()
        method_names = [row[0] for row in results]
        return method_names

    @open_and_close_sqlite3_connection
    def get_allele_names(self, method_name, connection=None):
        """ @brief: Returns a list of allele names given a \a method_name.
            @author: Dorjee
            @since: 2015-10-07
        """
        # Make method_name's passed case-insensitive.
        method_name = method_name.lower()

        if method_name not in self.get_method_names():
            raise ValueError('Invalid method_name: {}'.format(method_name))

        sql_where_clause = 'm.name = ?'
        sql_params = [method_name]
        sql = '''
            SELECT DISTINCT a.name
            FROM mhcnp_allele a JOIN mhcnp_method m ON a.method_id = m.id
            WHERE {0}
            ORDER BY a.name;
        '''.format(sql_where_clause)

        cur = connection.cursor()
        cur.execute(sql, sql_params)
        results = cur.fetchall()
        allele_names = [row[0] for row in results]
        return allele_names

    @open_and_close_sqlite3_connection
    def get_allowed_peptide_lengths(self, method_name, allele_name, connection=None):
        """ @brief: Returns a list of allele lengths given a \a method_name and a \a allele_name.
            @author: Dorjee
            @since: 2015-10-07
        """
        # Make method_name's passed case-insensitive.
        if method_name: method_name = method_name.lower()

        if method_name not in self.get_method_names():
            raise ValueError('Invalid method_name: {}'.format(method_name))
        if allele_name not in self.get_allele_names(method_name):
            raise ValueError('Invalid allele_name: {}'.format(allele_name))

        sql_where_clause = 'name = ?'
        sql_params = [allele_name]
        sql = '''
            SELECT DISTINCT length
            FROM mhcnp_allele
            WHERE {0}
            ORDER BY length;
        '''.format(sql_where_clause)

        cur = connection.cursor()
        cur.execute(sql, sql_params)
        results = cur.fetchall()
        lengths = [row[0] for row in results]
        return lengths


class NetCTLpanAlleleData(object):
    
    @open_and_close_sqlite3_connection
    def get_method_names(self, connection=None):
        """ @brief: Returns a list of valid method names for MHCNP.
            @author: Dorjee
            @since: 2016-08-05
        """
        sql = 'SELECT DISTINCT name FROM netctlpan_method ORDER BY name;'
        cur = connection.cursor()
        cur.execute(sql)
        results = cur.fetchall()
        method_names = [row[0] for row in results]
        return method_names
    
    @open_and_close_sqlite3_connection
    def get_species_list(self, connection=None):
        sql = 'SELECT DISTINCT species FROM netctlpan_allele ORDER BY species;'
        cur = connection.cursor()
        cur.execute(sql)
        results = cur.fetchall()
        species = [row[0] for row in results]
        return species
    
    @open_and_close_sqlite3_connection
    def get_allele_names_for_species(self, species, connection=None):
        """ @brief: Returns a list of allele names valid for use with the species indicated
                by \a species.
            @author: Dorjee
            @since: 
        """
        sql_where_clause = 'species = ?'
        sql_params = [species]
        sql = '''
            SELECT DISTINCT name FROM netctlpan_allele
            WHERE {0}
            ORDER BY name;
        '''.format(sql_where_clause)

        cur = connection.cursor()
        cur.execute(sql, sql_params)
        results = cur.fetchall()
        allele_names = [row[0] for row in results]
        return allele_names
    
    @open_and_close_sqlite3_connection
    def get_allowed_peptide_lengths(self, allele_name, connection=None):
        sql_where_clause = 'name = ?'
        sql_params = [allele_name]
        sql = '''
            SELECT DISTINCT length FROM netctlpan_allele
            WHERE {0}
            ORDER BY length;
        '''.format(sql_where_clause)

        cur = connection.cursor()
        cur.execute(sql, sql_params)
        results = cur.fetchall()
        lengths = [row[0] for row in results]
        return lengths
    
    @staticmethod
    def get_species_for_allele_name(allele_name):
        """ @brief: Returns \a species from \a given allele name.
            @author: Dorjee
            @since: 2016-11-01
        """
        # Remove all white-spaces
        allele_name = allele_name.strip()
        # Raise error if allele_name is not provided
        if not allele_name:
            raise ValueError('allele_name is a required keyword argument')

        species_indicator = allele_name[:2].lower()
        species_by_indicator = {
            'hl': 'human',
            'h-': 'mouse',
            'pa': 'chimpanzee',
            'ma': 'macaque',
            'go': 'gorilla',
            'sl': 'pig',
            'bo': 'cow',
            'rt': 'rat',
        }
        species = species_by_indicator[species_indicator]
        return species
    
def is_user_defined_allele(allele):
    ''' | *brief*: Returns True if *allele* is a user-defined sequence, False if it is an
        |    allele name.
        | *author*: Jivan
        | *created*: 2015-11-19
    '''
    if '-' in allele or '*' in allele or '1' in allele:
        ret = False
    else:
        ret = True

    return ret
