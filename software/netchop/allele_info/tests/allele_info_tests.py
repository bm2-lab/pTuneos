'''
Created on Sep 30, 2015

@author: jivan
'''
import os, json
from unittest import TestCase
import unittest

from allele_info import MHCIAlleleData, MHCIIAlleleData, MHCNPAlleleData, NetCTLpanAlleleData, is_user_defined_allele


class MHCIAlleleDataTests(TestCase):
    def get_allele_names_for_method_no_length_tests(self):
        miad = MHCIAlleleData()
        expected_allele_names_by_method = {
            'smm': ['BoLA-AW10', 'BoLA-D18.4', 'BoLA-HD6', 'BoLA-JSP.1', 'BoLA-T2C', 'BoLA-T2a', 'BoLA-T2b', 'H-2-Db', 'H-2-Dd', 'H-2-Kb', 'H-2-Kd', 'H-2-Kk', 'H-2-Ld', 'HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:02', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*02:11', 'HLA-A*02:12', 'HLA-A*02:16', 'HLA-A*02:17', 'HLA-A*02:19', 'HLA-A*02:50', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*24:03', 'HLA-A*25:01', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*29:02', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*32:07', 'HLA-A*32:15', 'HLA-A*33:01', 'HLA-A*66:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*68:23', 'HLA-A*69:01', 'HLA-A*80:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*14:02', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:09', 'HLA-B*15:17', 'HLA-B*15:42', 'HLA-B*18:01', 'HLA-B*27:05', 'HLA-B*27:20', 'HLA-B*35:01', 'HLA-B*35:03', 'HLA-B*38:01', 'HLA-B*39:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*40:13', 'HLA-B*42:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*45:01', 'HLA-B*46:01', 'HLA-B*48:01', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*54:01', 'HLA-B*57:01', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*73:01', 'HLA-B*83:01', 'HLA-C*03:03', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*08:02', 'HLA-C*12:03', 'HLA-C*14:02', 'HLA-C*15:02', 'HLA-E*01:01', 'HLA-E*01:03', 'Mamu-A*01', 'Mamu-A*02', 'Mamu-A*07', 'Mamu-A*11', 'Mamu-A*2201', 'Mamu-A*2601', 'Mamu-A2*0102', 'Mamu-A7*0103', 'Mamu-B*01', 'Mamu-B*03', 'Mamu-B*08', 'Mamu-B*1001', 'Mamu-B*17', 'Mamu-B*3901', 'Mamu-B*52', 'Mamu-B*6601', 'Mamu-B*8301', 'Mamu-B*8701', 'Patr-A*0101', 'Patr-A*0301', 'Patr-A*0401', 'Patr-A*0701', 'Patr-A*0901', 'Patr-B*0101', 'Patr-B*1301', 'Patr-B*2401', 'RT1A', 'SLA-1*0401', 'SLA-2*0401', 'SLA-3*0401'],
            'ann': ['BoLA-AW10', 'BoLA-D18.4', 'BoLA-HD6', 'BoLA-JSP.1', 'BoLA-T2C', 'BoLA-T2a', 'BoLA-T2b', 'H-2-Db', 'H-2-Dd', 'H-2-Kb', 'H-2-Kd', 'H-2-Kk', 'H-2-Ld', 'HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:02', 'HLA-A*02:03', 'HLA-A*02:05', 'HLA-A*02:06', 'HLA-A*02:07', 'HLA-A*02:11', 'HLA-A*02:12', 'HLA-A*02:16', 'HLA-A*02:17', 'HLA-A*02:19', 'HLA-A*02:50', 'HLA-A*03:01', 'HLA-A*03:02', 'HLA-A*03:19', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*24:03', 'HLA-A*25:01', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*29:02', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*32:07', 'HLA-A*32:15', 'HLA-A*33:01', 'HLA-A*66:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*68:23', 'HLA-A*69:01', 'HLA-A*80:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*14:01', 'HLA-B*14:02', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:09', 'HLA-B*15:17', 'HLA-B*18:01', 'HLA-B*27:05', 'HLA-B*27:20', 'HLA-B*35:01', 'HLA-B*35:03', 'HLA-B*37:01', 'HLA-B*38:01', 'HLA-B*39:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*40:13', 'HLA-B*42:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*45:01', 'HLA-B*45:06', 'HLA-B*46:01', 'HLA-B*48:01', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*54:01', 'HLA-B*57:01', 'HLA-B*57:03', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*73:01', 'HLA-B*81:01', 'HLA-B*83:01', 'HLA-C*03:03', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*08:02', 'HLA-C*12:03', 'HLA-C*14:02', 'HLA-C*15:02', 'HLA-E*01:01', 'HLA-E*01:03', 'Mamu-A*01', 'Mamu-A*02', 'Mamu-A*07', 'Mamu-A*11', 'Mamu-A*2201', 'Mamu-A*2601', 'Mamu-A2*0102', 'Mamu-A7*0103', 'Mamu-B*01', 'Mamu-B*03', 'Mamu-B*08', 'Mamu-B*1001', 'Mamu-B*17', 'Mamu-B*3901', 'Mamu-B*52', 'Mamu-B*6601', 'Mamu-B*8301', 'Mamu-B*8701', 'Patr-A*0101', 'Patr-A*0301', 'Patr-A*0401', 'Patr-A*0701', 'Patr-A*0901', 'Patr-B*0101', 'Patr-B*1301', 'Patr-B*2401', 'SLA-1*0401', 'SLA-1*0701', 'SLA-2*0401', 'SLA-3*0401'],
            'consensus': ['BoLA-AW10', 'BoLA-D18.4', 'BoLA-HD6', 'BoLA-JSP.1', 'BoLA-T2C', 'BoLA-T2a', 'BoLA-T2b', 'H-2-Db', 'H-2-Dd', 'H-2-Kb', 'H-2-Kd', 'H-2-Kk', 'H-2-Ld', 'HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:02', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*02:11', 'HLA-A*02:12', 'HLA-A*02:16', 'HLA-A*02:17', 'HLA-A*02:19', 'HLA-A*02:50', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*24:03', 'HLA-A*25:01', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*29:02', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*32:07', 'HLA-A*32:15', 'HLA-A*33:01', 'HLA-A*66:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*68:23', 'HLA-A*69:01', 'HLA-A*80:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*14:02', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:09', 'HLA-B*15:17', 'HLA-B*15:42', 'HLA-B*18:01', 'HLA-B*27:05', 'HLA-B*27:20', 'HLA-B*35:01', 'HLA-B*35:03', 'HLA-B*38:01', 'HLA-B*39:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*40:13', 'HLA-B*42:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*45:01', 'HLA-B*46:01', 'HLA-B*48:01', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*54:01', 'HLA-B*57:01', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*73:01', 'HLA-B*83:01', 'HLA-C*03:03', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*08:02', 'HLA-C*12:03', 'HLA-C*14:02', 'HLA-C*15:02', 'HLA-E*01:01', 'HLA-E*01:03', 'Mamu-A*01', 'Mamu-A*02', 'Mamu-A*07', 'Mamu-A*11', 'Mamu-A*2201', 'Mamu-A*2601', 'Mamu-A2*0102', 'Mamu-A7*0103', 'Mamu-B*01', 'Mamu-B*03', 'Mamu-B*08', 'Mamu-B*1001', 'Mamu-B*17', 'Mamu-B*3901', 'Mamu-B*52', 'Mamu-B*6601', 'Mamu-B*8301', 'Mamu-B*8701', 'Patr-A*0101', 'Patr-A*0301', 'Patr-A*0401', 'Patr-A*0701', 'Patr-A*0901', 'Patr-B*0101', 'Patr-B*1301', 'Patr-B*2401', 'RT1A', 'SLA-1*0401', 'SLA-2*0401', 'SLA-3*0401'],
        }

        for method_name, expected_allele_names in expected_allele_names_by_method.items():
            allele_names = miad.get_allele_names_for_method(method_name)
            self.assertEqual(allele_names, expected_allele_names)

        # method_name is a required argument
        self.assertRaises(TypeError, miad.get_allele_names_for_method)

    def get_allele_names_for_method_with_length_tests(self):
        miad = MHCIAlleleData()
        expected_allele_names_by_method_and_length = {
            ('smm', 8): [u'H-2-Db', u'H-2-Kb', u'H-2-Kd', u'H-2-Kk', u'HLA-A*01:01', u'HLA-A*02:01', u'HLA-A*02:02', u'HLA-A*02:03', u'HLA-A*02:06', u'HLA-A*03:01', u'HLA-A*11:01', u'HLA-A*23:01', u'HLA-A*24:02', u'HLA-A*26:01', u'HLA-A*29:02', u'HLA-A*30:02', u'HLA-A*68:02', u'HLA-B*07:02', u'HLA-B*08:01', u'HLA-B*18:01', u'HLA-B*27:05', u'HLA-B*35:01', u'HLA-B*40:01', u'HLA-B*40:02', u'HLA-B*44:02', u'HLA-B*44:03', u'HLA-B*45:01', u'HLA-B*51:01', u'HLA-B*53:01', u'HLA-B*54:01', u'Mamu-A*01', u'Mamu-A*02', u'Mamu-A*07', u'Mamu-A*11', u'Mamu-B*01', u'Mamu-B*03', u'Mamu-B*08', u'Mamu-B*17', u'Mamu-B*3901', u'Mamu-B*52'],
            ('ann', 10): ['BoLA-AW10', 'BoLA-D18.4', 'BoLA-HD6', 'BoLA-JSP.1', 'BoLA-T2C', 'BoLA-T2a', 'BoLA-T2b', 'H-2-Db', 'H-2-Dd', 'H-2-Kb', 'H-2-Kd', 'H-2-Kk', 'H-2-Ld', 'HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:02', 'HLA-A*02:03', 'HLA-A*02:05', 'HLA-A*02:06', 'HLA-A*02:07', 'HLA-A*02:11', 'HLA-A*02:12', 'HLA-A*02:16', 'HLA-A*02:17', 'HLA-A*02:19', 'HLA-A*02:50', 'HLA-A*03:01', 'HLA-A*03:02', 'HLA-A*03:19', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*24:03', 'HLA-A*25:01', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*29:02', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*32:07', 'HLA-A*32:15', 'HLA-A*33:01', 'HLA-A*66:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*68:23', 'HLA-A*69:01', 'HLA-A*80:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*14:01', 'HLA-B*14:02', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:09', 'HLA-B*15:17', 'HLA-B*18:01', 'HLA-B*27:05', 'HLA-B*27:20', 'HLA-B*35:01', 'HLA-B*35:03', 'HLA-B*37:01', 'HLA-B*38:01', 'HLA-B*39:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*40:13', 'HLA-B*42:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*45:01', 'HLA-B*45:06', 'HLA-B*46:01', 'HLA-B*48:01', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*54:01', 'HLA-B*57:01', 'HLA-B*57:03', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*73:01', 'HLA-B*81:01', 'HLA-B*83:01', 'HLA-C*03:03', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*08:02', 'HLA-C*12:03', 'HLA-C*14:02', 'HLA-C*15:02', 'HLA-E*01:01', 'HLA-E*01:03', 'Mamu-A*01', 'Mamu-A*02', 'Mamu-A*07', 'Mamu-A*11', 'Mamu-A*2201', 'Mamu-A*2601', 'Mamu-A2*0102', 'Mamu-A7*0103', 'Mamu-B*01', 'Mamu-B*03', 'Mamu-B*08', 'Mamu-B*1001', 'Mamu-B*17', 'Mamu-B*3901', 'Mamu-B*52', 'Mamu-B*6601', 'Mamu-B*8301', 'Mamu-B*8701', 'Patr-A*0101', 'Patr-A*0301', 'Patr-A*0401', 'Patr-A*0701', 'Patr-A*0901', 'Patr-B*0101', 'Patr-B*1301', 'Patr-B*2401', 'SLA-1*0401', 'SLA-1*0701', 'SLA-2*0401', 'SLA-3*0401'],
            ('consensus', 14): [u'BoLA-D18.4', u'BoLA-HD6', u'BoLA-JSP.1', u'BoLA-T2C', u'BoLA-T2a', u'BoLA-T2b', u'H-2-Db', u'H-2-Dd', u'H-2-Kb', u'H-2-Kd', u'H-2-Kk', u'H-2-Ld', u'HLA-A*01:01', u'HLA-A*02:01', u'HLA-A*02:02', u'HLA-A*02:03', u'HLA-A*02:06', u'HLA-A*02:11', u'HLA-A*02:12', u'HLA-A*02:16', u'HLA-A*02:17', u'HLA-A*02:19', u'HLA-A*02:50', u'HLA-A*03:01', u'HLA-A*11:01', u'HLA-A*23:01', u'HLA-A*24:02', u'HLA-A*24:03', u'HLA-A*25:01', u'HLA-A*26:01', u'HLA-A*26:02', u'HLA-A*26:03', u'HLA-A*29:02', u'HLA-A*30:01', u'HLA-A*30:02', u'HLA-A*31:01', u'HLA-A*32:01', u'HLA-A*32:07', u'HLA-A*32:15', u'HLA-A*33:01', u'HLA-A*66:01', u'HLA-A*68:01', u'HLA-A*68:02', u'HLA-A*68:23', u'HLA-A*69:01', u'HLA-A*80:01', u'HLA-B*07:02', u'HLA-B*08:01', u'HLA-B*08:02', u'HLA-B*08:03', u'HLA-B*14:02', u'HLA-B*15:01', u'HLA-B*15:02', u'HLA-B*15:03', u'HLA-B*15:09', u'HLA-B*15:17', u'HLA-B*18:01', u'HLA-B*27:05', u'HLA-B*27:20', u'HLA-B*35:01', u'HLA-B*35:03', u'HLA-B*38:01', u'HLA-B*39:01', u'HLA-B*40:01', u'HLA-B*40:02', u'HLA-B*40:13', u'HLA-B*42:01', u'HLA-B*44:02', u'HLA-B*44:03', u'HLA-B*45:01', u'HLA-B*46:01', u'HLA-B*48:01', u'HLA-B*51:01', u'HLA-B*53:01', u'HLA-B*54:01', u'HLA-B*57:01', u'HLA-B*58:01', u'HLA-B*73:01', u'HLA-B*83:01', u'HLA-C*03:03', u'HLA-C*04:01', u'HLA-C*05:01', u'HLA-C*06:02', u'HLA-C*07:01', u'HLA-C*07:02', u'HLA-C*08:02', u'HLA-C*12:03', u'HLA-C*14:02', u'HLA-C*15:02', u'HLA-E*01:01', u'Mamu-A*01', u'Mamu-A*02', u'Mamu-A*07', u'Mamu-A*11', u'Mamu-A*2201', u'Mamu-A*2601', u'Mamu-A2*0102', u'Mamu-A7*0103', u'Mamu-B*01', u'Mamu-B*03', u'Mamu-B*08', u'Mamu-B*1001', u'Mamu-B*17', u'Mamu-B*3901', u'Mamu-B*52', u'Mamu-B*6601', u'Mamu-B*8301', u'Mamu-B*8701', u'Patr-A*0101', u'Patr-A*0301', u'Patr-A*0401', u'Patr-A*0701', u'Patr-A*0901', u'Patr-B*0101', u'Patr-B*1301', u'Patr-B*2401', u'SLA-1*0401', u'SLA-2*0401', u'SLA-3*0401'],
        }

        for inargs, expected_allele_names in expected_allele_names_by_method_and_length.items():
            method_name, binding_length = inargs
            allele_names = miad.get_allele_names_for_method(
                                method_name, binding_length=binding_length)
            msg = "({}, {}): {} != {}".format(
                    method_name, binding_length, allele_names, expected_allele_names)
            self.assertEqual(allele_names, expected_allele_names, msg=msg)

        # method_name is a required argument
        self.assertRaises(TypeError, miad.get_allele_names_for_method)

    def get_all_allele_names_tests(self):
        miad = MHCIAlleleData()
        expected_names_dir = os.path.dirname(os.path.realpath(__file__))
        expected_names_filepath = os.path.join(expected_names_dir, 'all_allele_names.json')
        with open(expected_names_filepath) as f:
            expected_allele_names = json.loads(f.read())
        allele_names = miad.get_all_allele_names()
        print('Resulting allele names:\n{}'.format(allele_names))
        self.assertEqual(allele_names, expected_allele_names)

    def get_allele_names_tests(self):
        miad = MHCIAlleleData()
        # Filter only by species
        species = 'chimpanzee'
        expected_allele_names = ['Patr-A*0101', 'Patr-A*0201', 'Patr-A*0301', 'Patr-A*0302', 'Patr-A*0401', 'Patr-A*0402', 'Patr-A*0404', 'Patr-A*0501', 'Patr-A*0601', 'Patr-A*0602', 'Patr-A*0701', 'Patr-A*0801', 'Patr-A*0802', 'Patr-A*0803', 'Patr-A*0901', 'Patr-A*0902', 'Patr-A*1001', 'Patr-A*1101', 'Patr-A*1201', 'Patr-A*1301', 'Patr-A*1401', 'Patr-A*1501', 'Patr-A*1502', 'Patr-A*1601', 'Patr-A*1701', 'Patr-A*1702', 'Patr-A*1703', 'Patr-A*1801', 'Patr-A*2301', 'Patr-A*2401', 'Patr-B*0101', 'Patr-B*0102', 'Patr-B*0201', 'Patr-B*0203', 'Patr-B*0301', 'Patr-B*0302', 'Patr-B*0401', 'Patr-B*0402', 'Patr-B*0501', 'Patr-B*0502', 'Patr-B*0601', 'Patr-B*0701', 'Patr-B*0702', 'Patr-B*0801', 'Patr-B*0802', 'Patr-B*0901', 'Patr-B*1001', 'Patr-B*1101', 'Patr-B*1102', 'Patr-B*1202', 'Patr-B*1301', 'Patr-B*1401', 'Patr-B*1601', 'Patr-B*1602', 'Patr-B*1701', 'Patr-B*1702', 'Patr-B*1703', 'Patr-B*1801', 'Patr-B*1901', 'Patr-B*2001', 'Patr-B*2101', 'Patr-B*2201', 'Patr-B*2202', 'Patr-B*2301', 'Patr-B*2302', 'Patr-B*2303', 'Patr-B*2401', 'Patr-B*2402', 'Patr-B*2501', 'Patr-B*2601', 'Patr-B*2701', 'Patr-B*2801', 'Patr-B*2901', 'Patr-B*3001', 'Patr-B*3501', 'Patr-B*3601', 'Patr-B*3701', 'Patr-C*0201', 'Patr-C*0202', 'Patr-C*0203', 'Patr-C*0204', 'Patr-C*0205', 'Patr-C*0206', 'Patr-C*0301', 'Patr-C*0302', 'Patr-C*0303', 'Patr-C*0304', 'Patr-C*0401', 'Patr-C*0501', 'Patr-C*0502', 'Patr-C*0601', 'Patr-C*0701', 'Patr-C*0801', 'Patr-C*0901', 'Patr-C*0902', 'Patr-C*0903', 'Patr-C*0904', 'Patr-C*0905', 'Patr-C*1001', 'Patr-C*1101', 'Patr-C*1201', 'Patr-C*1301', 'Patr-C*1302', 'Patr-C*1501', 'Patr-C*1601', ]
        allele_names = miad.get_allele_names(species=species)
        self.assertEqual(allele_names, expected_allele_names)

        # Filter by species & method
        species = 'chimpanzee'
        method_name = 'netmhcpan'
        expected_allele_names = ['Patr-A*0101', 'Patr-A*0201', 'Patr-A*0301', 'Patr-A*0302', 'Patr-A*0401', 'Patr-A*0402', 'Patr-A*0404', 'Patr-A*0501', 'Patr-A*0601', 'Patr-A*0602', 'Patr-A*0701', 'Patr-A*0801', 'Patr-A*0802', 'Patr-A*0803', 'Patr-A*0901', 'Patr-A*0902', 'Patr-A*1001', 'Patr-A*1101', 'Patr-A*1201', 'Patr-A*1301', 'Patr-A*1401', 'Patr-A*1501', 'Patr-A*1502', 'Patr-A*1601', 'Patr-A*1701', 'Patr-A*1702', 'Patr-A*1703', 'Patr-A*1801', 'Patr-A*2301', 'Patr-A*2401', 'Patr-B*0101', 'Patr-B*0102', 'Patr-B*0201', 'Patr-B*0203', 'Patr-B*0301', 'Patr-B*0302', 'Patr-B*0401', 'Patr-B*0402', 'Patr-B*0501', 'Patr-B*0502', 'Patr-B*0601', 'Patr-B*0701', 'Patr-B*0702', 'Patr-B*0801', 'Patr-B*0802', 'Patr-B*0901', 'Patr-B*1001', 'Patr-B*1101', 'Patr-B*1102', 'Patr-B*1202', 'Patr-B*1301', 'Patr-B*1401', 'Patr-B*1601', 'Patr-B*1602', 'Patr-B*1701', 'Patr-B*1702', 'Patr-B*1703', 'Patr-B*1801', 'Patr-B*1901', 'Patr-B*2001', 'Patr-B*2101', 'Patr-B*2201', 'Patr-B*2202', 'Patr-B*2301', 'Patr-B*2302', 'Patr-B*2303', 'Patr-B*2401', 'Patr-B*2402', 'Patr-B*2501', 'Patr-B*2601', 'Patr-B*2701', 'Patr-B*2801', 'Patr-B*2901', 'Patr-B*3001', 'Patr-B*3501', 'Patr-B*3601', 'Patr-B*3701', 'Patr-C*0201', 'Patr-C*0202', 'Patr-C*0203', 'Patr-C*0204', 'Patr-C*0205', 'Patr-C*0206', 'Patr-C*0301', 'Patr-C*0302', 'Patr-C*0303', 'Patr-C*0304', 'Patr-C*0401', 'Patr-C*0501', 'Patr-C*0502', 'Patr-C*0601', 'Patr-C*0701', 'Patr-C*0801', 'Patr-C*0901', 'Patr-C*0902', 'Patr-C*0903', 'Patr-C*0904', 'Patr-C*0905', 'Patr-C*1001', 'Patr-C*1101', 'Patr-C*1201', 'Patr-C*1301', 'Patr-C*1302', 'Patr-C*1501', 'Patr-C*1601', ]
        allele_names = miad.get_allele_names(species=species, method_name=method_name)
        self.assertEqual(allele_names, expected_allele_names)

        # Filter by species & frequency - checking frequency_cutoff=1, species=human
        species = 'human'
        frequency_cutoff = 1.0
        expected_allele_names = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*25:01', 'HLA-A*26:01', 'HLA-A*29:02', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*33:01', 'HLA-A*33:03', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*74:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*13:01', 'HLA-B*13:02', 'HLA-B*14:02', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:25', 'HLA-B*18:01', 'HLA-B*27:02', 'HLA-B*27:05', 'HLA-B*35:01', 'HLA-B*35:03', 'HLA-B*37:01', 'HLA-B*38:01', 'HLA-B*39:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*46:01', 'HLA-B*48:01', 'HLA-B*49:01', 'HLA-B*50:01', 'HLA-B*51:01', 'HLA-B*52:01', 'HLA-B*53:01', 'HLA-B*55:01', 'HLA-B*56:01', 'HLA-B*57:01', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-C*01:02', 'HLA-C*02:02', 'HLA-C*02:09', 'HLA-C*03:02', 'HLA-C*03:03', 'HLA-C*03:04', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*07:04', 'HLA-C*08:01', 'HLA-C*08:02', 'HLA-C*12:02', 'HLA-C*12:03', 'HLA-C*14:02', 'HLA-C*15:02', 'HLA-C*16:01', 'HLA-C*17:01', 'HLA-E*01:01', 'HLA-E*01:03', 'HLA-G*01:01', 'HLA-G*01:02', 'HLA-G*01:03', 'HLA-G*01:04', 'HLA-G*01:06']
        allele_names = miad.get_allele_names(
                          species=species, frequency_cutoff=frequency_cutoff)
        self.assertEqual(allele_names, expected_allele_names)

        # Call without any parameters (species is required)
        self.assertRaises(ValueError, miad.get_allele_names, *[miad])

        # Filter by species, method, & frequency
        species = 'human'
        method_name = 'smm'
        frequency_cutoff = 4.0
        expected_allele_names = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*24:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*18:01', 'HLA-B*35:01', 'HLA-B*40:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*51:01', 'HLA-C*03:03', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*12:03', 'HLA-E*01:01']
        allele_names = miad.get_allele_names(
                          species=species, method_name=method_name, frequency_cutoff=frequency_cutoff)
        self.assertEqual(allele_names, expected_allele_names)

        # Call without any parameters (species is required)
        self.assertRaises(ValueError, miad.get_allele_names, *[miad])

        # Call only with method (species is required)
        self.assertRaises(ValueError, miad.get_allele_names, *[], **{'method_name': method_name})

    def get_allele_frequency_tests(self):
        miad = MHCIAlleleData()
        name_frequency = [('HLA-A*01:17', None),
                          ('HLA-A*01:01', 15.9), ('HLA-A*01:06', 0.01),
                          ('HLA-B*15:27', 0.28), ('HLA-C*08:04', 0.41), ]
        names = [t[0] for t in name_frequency]
        expected_frequencies = [t[1] for t in name_frequency]
        frequencies = miad.get_allele_frequencies(names)
        self.assertEqual(frequencies, expected_frequencies)

    def get_method_names_tests(self):
        miad = MHCIAlleleData()
        # Some methods are hidden from users but still available via an optional keyword.
        expected_method_names_with_hidden = [
            'ann', 'arb', 'comblib_sidney2008', 'consensus', 'netmhccons',
            'netmhcpan', 'pickpocket', 'recommended', 'smm', 'smmpmbec',
        ]
        expected_method_names_without_hidden = [
            'ann', 'comblib_sidney2008', 'consensus', 'netmhccons',
            'netmhcpan', 'pickpocket', 'recommended', 'smm', 'smmpmbec',
        ]

        # Call without any filter parameters
        method_names_without_hidden = miad.get_method_names()
        self.assertEqual(method_names_without_hidden, expected_method_names_without_hidden)
        method_names_with_hidden = miad.get_method_names(include_hidden=True)
        self.assertEqual(method_names_with_hidden, expected_method_names_with_hidden)

        # Call with allele_name and binding_length parameters ()
        allele_name = 'HLA-A*02:01'
        binding_length = 10
        method_names = miad.get_method_names(allele_name=allele_name, binding_length=binding_length)
        for method_name in method_names:
            self.assertIn(method_name, expected_method_names_without_hidden)

    def get_species_list_tests(self):
        miad = MHCIAlleleData()
        expected_species_by_method = {
            'ann': ['chimpanzee', 'cow', 'human', 'macaque', 'mouse', 'pig'],
            'comblib_sidney2008': ['human', 'mouse'],
            'consensus': ['chimpanzee', 'cow', 'human', 'macaque', 'mouse', 'pig', 'rat'],
            'netmhccons': ['chimpanzee', 'cow', 'gorilla', 'human', 'macaque', 'mouse', 'pig'],
            'netmhcpan': ['chimpanzee', 'cow', 'gorilla', 'human', 'macaque', 'mouse', 'pig'],
            'pickpocket': ['chimpanzee', 'cow', 'gorilla', 'human', 'macaque', 'mouse', 'pig'],
            'recommended': ['chimpanzee', 'cow', 'gorilla', 'human', 'macaque', 'mouse', 'pig', 'rat'],
            'smm': ['chimpanzee', 'cow', 'human', 'macaque', 'mouse', 'pig', 'rat'],
            'smmpmbec': ['chimpanzee', 'cow', 'human', 'macaque', 'mouse', 'pig', 'rat'],
        }
        for method_name, expected_species in expected_species_by_method.items():
            species = miad.get_species_list(method_name=method_name)
            self.assertEqual(species, expected_species)

        expected_species = ['chimpanzee', 'cow', 'gorilla', 'human', 'macaque', 'mouse', 'pig', 'rat']
        species_list = miad.get_species_list()
        self.assertEqual(species_list, expected_species)

    def get_species_for_allele_name_tests(self):
        allele_name = 'HLA-A*02:01'
        miad = MHCIAlleleData()
        species = miad.get_species_for_allele_name(allele_name)
        self.assertEqual(species, 'human')

        allele_name = 'Patr-A*0101'
        miad = MHCIAlleleData()
        species = miad.get_species_for_allele_name(allele_name)
        self.assertEqual(species, 'chimpanzee')

    def get_allowed_peptide_lengths_tests(self):
        miad = MHCIAlleleData()
        expected_lengths_by_method_and_allele = {
            ('recommended', 'SLA-6*0104'): [8, 9, 10, 11, 12, 13, 14],
            ('consensus', 'Patr-A*0901'): [8, 9, 10, 11, 12, 13, 14],
            ('netmhcpan', 'BoLA-AW10'): [8, 9, 10, 11, 12, 13, 14],
            ('ann', 'BoLA-JSP.1'): [8, 9, 10, 11, 12, 13, 14],
            ('smmpmbec', 'H-2-Db'): [8, 9, 10, 11],
            ('smm', 'HLA-A*02:01'): [8, 9, 10, 11],
        }

        for (method_name, allele_name), expected_lengths in expected_lengths_by_method_and_allele.items():
            lengths = miad.get_allowed_peptide_lengths(method_name=method_name, allele_name=allele_name)
            msg = '{} != {} for {}:{}'.format(lengths, expected_lengths, method_name, allele_name)
            self.assertEqual(lengths, expected_lengths, msg)

    def get_reference_set_tests(self):
        miad = MHCIAlleleData()
        expected_reference_set = [('HLA-A*01:01', '9'), ('HLA-A*01:01', '10'), ('HLA-A*02:01', '9'), ('HLA-A*02:01', '10'), ('HLA-A*02:03', '9'), ('HLA-A*02:03', '10'), ('HLA-A*02:06', '9'), ('HLA-A*02:06', '10'), ('HLA-A*03:01', '9'), ('HLA-A*03:01', '10'), ('HLA-A*11:01', '9'), ('HLA-A*11:01', '10'), ('HLA-A*23:01', '9'), ('HLA-A*23:01', '10'), ('HLA-A*24:02', '9'), ('HLA-A*24:02', '10'), ('HLA-A*26:01', '9'), ('HLA-A*26:01', '10'), ('HLA-A*30:01', '9'), ('HLA-A*30:01', '10'), ('HLA-A*30:02', '9'), ('HLA-A*30:02', '10'), ('HLA-A*31:01', '9'), ('HLA-A*31:01', '10'), ('HLA-A*32:01', '9'), ('HLA-A*32:01', '10'), ('HLA-A*33:01', '9'), ('HLA-A*33:01', '10'), ('HLA-A*68:01', '9'), ('HLA-A*68:01', '10'), ('HLA-A*68:02', '9'), ('HLA-A*68:02', '10'), ('HLA-B*07:02', '9'), ('HLA-B*07:02', '10'), ('HLA-B*08:01', '9'), ('HLA-B*08:01', '10'), ('HLA-B*15:01', '9'), ('HLA-B*15:01', '10'), ('HLA-B*35:01', '9'), ('HLA-B*35:01', '10'), ('HLA-B*40:01', '9'), ('HLA-B*40:01', '10'), ('HLA-B*44:02', '9'), ('HLA-B*44:02', '10'), ('HLA-B*44:03', '9'), ('HLA-B*44:03', '10'), ('HLA-B*51:01', '9'), ('HLA-B*51:01', '10'), ('HLA-B*53:01', '9'), ('HLA-B*53:01', '10'), ('HLA-B*57:01', '9'), ('HLA-B*57:01', '10'), ('HLA-B*58:01', '9'), ('HLA-B*58:01', '10')]
        reference_set = miad.get_reference_set()
        self.assertEqual(reference_set, expected_reference_set)


class MHCIIAlleleDataTests(TestCase):

    def get_allele_names_tests(self):
        miiad = MHCIIAlleleData()
        # Filter by method and locus
        method_name = 'smm_align'
        locus_name = 'DR'
        expected_DR_allele_names = ["DRB1*01:01", "DRB1*03:01", "DRB1*04:01", "DRB1*04:04", "DRB1*04:05", "DRB1*07:01", "DRB1*08:02", "DRB1*09:01", "DRB1*11:01", "DRB1*12:01", "DRB1*13:02", "DRB1*15:01", "DRB3*01:01", "DRB4*01:01", "DRB5*01:01"]
        DR_allele_names = miiad.get_allele_names(method_name=method_name, locus_name=locus_name)
        self.assertEqual(DR_allele_names, expected_DR_allele_names)

        expected_smm_align_allele_names = ['DPA1*01/DPB1*04:01', 'DPA1*01:03/DPB1*02:01', 'DPA1*02:01/DPB1*01:01', 'DPA1*02:01/DPB1*05:01', 'DPA1*03:01/DPB1*04:02', 'DQA1*01:01/DQB1*05:01', 'DQA1*01:02/DQB1*06:02', 'DQA1*03:01/DQB1*03:02', 'DQA1*04:01/DQB1*04:02', 'DQA1*05:01/DQB1*02:01', 'DQA1*05:01/DQB1*03:01', 'DRB1*01:01', 'DRB1*03:01', 'DRB1*04:01', 'DRB1*04:04', 'DRB1*04:05', 'DRB1*07:01', 'DRB1*08:02', 'DRB1*09:01', 'DRB1*11:01', 'DRB1*12:01', 'DRB1*13:02', 'DRB1*15:01', 'DRB3*01:01', 'DRB4*01:01', 'DRB5*01:01', 'H2-IAb', 'H2-IAd', 'H2-IEd']
        smm_align_allele_names = miiad.get_allele_names(method_name=method_name)
        self.assertEqual(smm_align_allele_names, expected_smm_align_allele_names)
        
        # Call without any parameters (method and locus are required)
        self.assertRaises(ValueError, miiad.get_allele_names, *[])

    def get_method_names_tests(self):
        miiad = MHCIIAlleleData()
        expected_method_names = ["comblib", "consensus", "netmhciipan", "nn_align", "recommended", "smm_align", "tepitope"]
        method_names = miiad.get_method_names()
        self.assertEqual(method_names, expected_method_names)

    def get_locus_names_tests(self):
        miiad = MHCIIAlleleData()
        expected_locus_names = ["DP", "DQ", "DR", "H2"]
        locus_names = miiad.get_locus_names()
        self.assertEqual(locus_names, expected_locus_names)

    def get_locus_names_for_method_tests(self):    
        miiad = MHCIIAlleleData()
        # Filter by method
        method_name = 'comblib'
        expected_locus_names = ["DP", "DQ", "DR"]
        locus_names = miiad.get_locus_names_for_method(method_name=method_name)
        self.assertEqual(locus_names, expected_locus_names)
    
    def get_alpha_chain_tests(self):
        miiad = MHCIIAlleleData()
        # Filter by method and allele
        method_name = 'recommended'
        locus_name = 'DP'
        expected_alpha_chains = ['DPA1*01', 'DPA1*01:03', 'DPA1*01:04', 'DPA1*01:05', 'DPA1*01:06', 'DPA1*01:07', 'DPA1*01:08', 'DPA1*01:09', 'DPA1*01:10', 'DPA1*02:01', u'DPA1*02:02', 'DPA1*02:03', 'DPA1*02:04', 'DPA1*03:01', 'DPA1*03:02', 'DPA1*03:03', 'DPA1*04:01']
        alpha_chains = miiad.get_alpha_chain(method_name=method_name, locus_name=locus_name)
        self.assertEqual(alpha_chains, expected_alpha_chains)

        # Call without any parameters (method and allele are required)
        self.assertRaises(ValueError, miiad.get_beta_chain, *[])
            
    def get_beta_chain_tests(self):
        miiad = MHCIIAlleleData()
        # Filter by method and allele
        method_name = 'recommended'
        allele_name = 'DPA1*01:03'
        expected_beta_chains = ["DPB1*01:01", "DPB1*02:01", "DPB1*02:02", "DPB1*03:01", "DPB1*04:01", "DPB1*04:02", "DPB1*05:01", "DPB1*06:01", "DPB1*08:01", "DPB1*09:01", "DPB1*100:01", "DPB1*101:01", "DPB1*102:01", "DPB1*103:01", "DPB1*104:01", "DPB1*105:01", "DPB1*106:01", "DPB1*107:01", "DPB1*108:01", "DPB1*109:01", "DPB1*10:01", "DPB1*110:01", "DPB1*111:01", "DPB1*112:01", "DPB1*113:01", "DPB1*114:01", "DPB1*115:01", "DPB1*116:01", "DPB1*117:01", "DPB1*118:01", "DPB1*119:01", "DPB1*11:01", "DPB1*121:01", "DPB1*122:01", "DPB1*123:01", "DPB1*124:01", "DPB1*125:01", "DPB1*126:01", "DPB1*127:01", "DPB1*128:01", "DPB1*129:01", "DPB1*130:01", "DPB1*131:01", "DPB1*132:01", "DPB1*133:01", "DPB1*134:01", "DPB1*13:01", "DPB1*14:01", "DPB1*15:01", "DPB1*16:01", "DPB1*17:01", "DPB1*18:01", "DPB1*19:01", "DPB1*20:01", "DPB1*21:01", "DPB1*22:01", "DPB1*23:01", "DPB1*24:01", "DPB1*25:01", "DPB1*26:01", "DPB1*27:01", "DPB1*28:01", "DPB1*29:01", "DPB1*30:01", "DPB1*31:01", "DPB1*32:01", "DPB1*33:01", "DPB1*34:01", "DPB1*35:01", "DPB1*36:01", "DPB1*37:01", "DPB1*38:01", "DPB1*39:01", "DPB1*40:01", "DPB1*41:01", "DPB1*44:01", "DPB1*45:01", "DPB1*46:01", "DPB1*47:01", "DPB1*48:01", "DPB1*49:01", "DPB1*50:01", "DPB1*51:01", "DPB1*52:01", "DPB1*53:01", "DPB1*54:01", "DPB1*55:01", "DPB1*56:01", "DPB1*58:01", "DPB1*59:01", "DPB1*60:01", "DPB1*62:01", "DPB1*63:01", "DPB1*65:01", "DPB1*66:01", "DPB1*67:01", "DPB1*68:01", "DPB1*69:01", "DPB1*70:01", "DPB1*71:01", "DPB1*72:01", "DPB1*73:01", "DPB1*74:01", "DPB1*75:01", "DPB1*76:01", "DPB1*77:01", "DPB1*78:01", "DPB1*79:01", "DPB1*80:01", "DPB1*81:01", "DPB1*82:01", "DPB1*83:01", "DPB1*84:01", "DPB1*85:01", "DPB1*86:01", "DPB1*87:01", "DPB1*88:01", "DPB1*89:01", "DPB1*90:01", "DPB1*91:01", "DPB1*92:01", "DPB1*93:01", "DPB1*94:01", "DPB1*95:01", "DPB1*96:01", "DPB1*97:01", "DPB1*98:01", "DPB1*99:01"]
        beta_chains = miiad.get_beta_chain(method_name=method_name, allele_name=allele_name)
        self.assertEqual(beta_chains, expected_beta_chains)
        
        # Call without any parameters (method and allele are required)
        self.assertRaises(ValueError, miiad.get_beta_chain, *[])

    def get_alpha_chains_for_locus_tests(self):
        miiad = MHCIIAlleleData()
        # Filter by locus
        locus_name = 'DP'
        expected_alpha_chains = ["DPA1*01", "DPA1*01:03", "DPA1*01:04", "DPA1*01:05", "DPA1*01:06", "DPA1*01:07", "DPA1*01:08", "DPA1*01:09", "DPA1*01:10", "DPA1*02:01", "DPA1*02:02", "DPA1*02:03", "DPA1*02:04", "DPA1*03:01", "DPA1*03:02", "DPA1*03:03", "DPA1*04:01"] 
        alpha_chains = miiad.get_alpha_chains_for_locus(locus_name=locus_name)
        self.assertEqual(alpha_chains, expected_alpha_chains)

    def get_beta_chains_for_locus_tests(self):
        miiad = MHCIIAlleleData()
        # Filter by locus
        locus_name = 'DQ'
        expected_beta_chains = ["DQB1*02:01", "DQB1*02:02", "DQB1*02:03", "DQB1*02:04", "DQB1*02:05", "DQB1*02:06", "DQB1*03:01", "DQB1*03:02", "DQB1*03:03", "DQB1*03:04", "DQB1*03:05", "DQB1*03:06", "DQB1*03:07", "DQB1*03:08", "DQB1*03:09", "DQB1*03:10", "DQB1*03:11", "DQB1*03:12", "DQB1*03:13", "DQB1*03:14", "DQB1*03:15", "DQB1*03:16", "DQB1*03:17", "DQB1*03:18", "DQB1*03:19", "DQB1*03:20", "DQB1*03:21", "DQB1*03:22", "DQB1*03:23", "DQB1*03:24", "DQB1*03:25", "DQB1*03:26", "DQB1*03:27", "DQB1*03:28", "DQB1*03:29", "DQB1*03:30", "DQB1*03:31", "DQB1*03:32", "DQB1*03:33", "DQB1*03:34", "DQB1*03:35", "DQB1*03:36", "DQB1*03:37", "DQB1*03:38", "DQB1*04:01", "DQB1*04:02", "DQB1*04:03", "DQB1*04:04", "DQB1*04:05", "DQB1*04:06", "DQB1*04:07", "DQB1*04:08", "DQB1*05:01", "DQB1*05:02", "DQB1*05:03", "DQB1*05:05", "DQB1*05:06", "DQB1*05:07", "DQB1*05:08", "DQB1*05:09", "DQB1*05:10", "DQB1*05:11", "DQB1*05:12", "DQB1*05:13", "DQB1*05:14", "DQB1*06:01", "DQB1*06:02", "DQB1*06:03", "DQB1*06:04", "DQB1*06:07", "DQB1*06:08", "DQB1*06:09", "DQB1*06:10", "DQB1*06:11", "DQB1*06:12", "DQB1*06:14", "DQB1*06:15", "DQB1*06:16", "DQB1*06:17", "DQB1*06:18", "DQB1*06:19", "DQB1*06:21", "DQB1*06:22", "DQB1*06:23", "DQB1*06:24", "DQB1*06:25", "DQB1*06:27", "DQB1*06:28", "DQB1*06:29", "DQB1*06:30", "DQB1*06:31", "DQB1*06:32", "DQB1*06:33", "DQB1*06:34", "DQB1*06:35", "DQB1*06:36", "DQB1*06:37", "DQB1*06:38", "DQB1*06:39", "DQB1*06:40", "DQB1*06:41", "DQB1*06:42", "DQB1*06:43", "DQB1*06:44"]
        beta_chains = miiad.get_beta_chains_for_locus(locus_name=locus_name)
        self.assertEqual(beta_chains, expected_beta_chains)
    
    def get_reference_set_tests(self):
        miiad = MHCIIAlleleData()
        expected_reference_set = ["HLA-DPA1*01/DPB1*04:01", "HLA-DPA1*01:03/DPB1*02:01", "HLA-DPA1*02:01/DPB1*01:01", "HLA-DPA1*02:01/DPB1*05:01", "HLA-DPA1*03:01/DPB1*04:02", "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*04:01/DQB1*04:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*05:01/DQB1*03:01", "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*04:05", "HLA-DRB1*07:01", "HLA-DRB1*08:02", "HLA-DRB1*09:01", "HLA-DRB1*11:01", "HLA-DRB1*12:01", "HLA-DRB1*13:02", "HLA-DRB1*15:01", "HLA-DRB3*01:01", "HLA-DRB3*02:02", "HLA-DRB4*01:01", "HLA-DRB5*01:01"]
        reference_set = miiad.get_reference_set()
        self.assertEqual(reference_set, expected_reference_set)


class MHCNPAlleleDataTests(TestCase):

    def get_method_names_tests(self):
        npad = MHCNPAlleleData()
        expected_method_names = ['mhcnp']
        method_names = npad.get_method_names()
        self.assertEqual(method_names, expected_method_names)

    def get_allele_names_tests(self):
        npad = MHCNPAlleleData()
        # Filter by method
        method_name = 'mhcnp'
        expected_allele_names = ["H-2-Db", "H-2-Kb", "HLA-A*02:01", "HLA-B*07:02", "HLA-B*35:01", "HLA-B*44:03", "HLA-B*53:01", "HLA-B*57:01"]
        mhcnp_allele_names = npad.get_allele_names(method_name=method_name)
        self.assertEqual(mhcnp_allele_names, expected_allele_names)

    def get_allowed_peptide_lengths_tests(self):
        npad = MHCNPAlleleData()
        # Filter by method & allele
        method_name = 'mhcnp'
        allele_name = 'HLA-A*02:01'
        expected_peptide_lengths = [8, 9, 10, 11]
        peptide_lengths_0201 = npad.get_allowed_peptide_lengths(method_name=method_name, allele_name=allele_name)
        self.assertEqual(peptide_lengths_0201, expected_peptide_lengths)


class NetCTLpanAlleleDataTests(TestCase):
 
    def get_method_names_tests(self):
        ctlpanad = NetCTLpanAlleleData()
        expected_method_name = ['netctlpan']
        method_name = ctlpanad.get_method_names()
        self.assertEqual(method_name, expected_method_name)
    
    def get_species_list_tests(self):
        ctlpanad = NetCTLpanAlleleData()
        expected_species = ['chimpanzee', 'gorilla', 'human', 'macaque', 'mouse', 'pig']
        species_list = ctlpanad.get_species_list()
        self.assertEqual(species_list, expected_species)
        
    def get_allele_names_for_species_test(self):
        ctlpanad = NetCTLpanAlleleData()
        
        # Filter only by species
        species = 'chimpanzee'
        expected_allele_names = ['Patr-A0101', 'Patr-A0201', 'Patr-A0301', 'Patr-A0302', 'Patr-A0401', 'Patr-A0402', 'Patr-A0404', 'Patr-A0501', 'Patr-A0601', 'Patr-A0602', 'Patr-A0701', 'Patr-A0801', 'Patr-A0802', 'Patr-A0803', 'Patr-A0901', 'Patr-A0902', 'Patr-A1001', 'Patr-A1101', 'Patr-A1201', 'Patr-A1301', 'Patr-A1401', 'Patr-A1501', 'Patr-A1502', 'Patr-A1601', 'Patr-A1701', 'Patr-A1702', 'Patr-A1703', 'Patr-A1801', 'Patr-A2301', 'Patr-A2401', 'Patr-B0101', 'Patr-B0102', 'Patr-B0201', 'Patr-B0203', 'Patr-B0301', 'Patr-B0302', 'Patr-B0401', 'Patr-B0402', 'Patr-B0501', 'Patr-B0502', 'Patr-B0601', 'Patr-B0701', 'Patr-B0702', 'Patr-B0801', 'Patr-B0802', 'Patr-B0901', 'Patr-B1001', 'Patr-B1101', 'Patr-B1102', 'Patr-B1202', 'Patr-B1301', 'Patr-B1401', 'Patr-B1601', 'Patr-B1602', 'Patr-B1701', 'Patr-B1702', 'Patr-B1703', 'Patr-B1801', 'Patr-B1901', 'Patr-B2001', 'Patr-B2101', 'Patr-B2201', 'Patr-B2202', 'Patr-B2301', 'Patr-B2302', 'Patr-B2303', 'Patr-B2401', 'Patr-B2402', 'Patr-B2501', 'Patr-B2601', 'Patr-B2701', 'Patr-B2801', 'Patr-B2901', 'Patr-B3001', 'Patr-B3501', 'Patr-B3601', 'Patr-B3701', 'Patr-C0201', 'Patr-C0202', 'Patr-C0203', 'Patr-C0204', 'Patr-C0205', 'Patr-C0206', 'Patr-C0301', 'Patr-C0302', 'Patr-C0303', 'Patr-C0304', 'Patr-C0401', 'Patr-C0501', 'Patr-C0502', 'Patr-C0601', 'Patr-C0701', 'Patr-C0801', 'Patr-C0901', 'Patr-C0902', 'Patr-C0903', 'Patr-C0904', 'Patr-C0905', 'Patr-C1001', 'Patr-C1101', 'Patr-C1201', 'Patr-C1301', 'Patr-C1302', 'Patr-C1501', 'Patr-C1601']
        allele_names = ctlpanad.get_allele_names_for_species(species=species)
        self.assertEqual(allele_names, expected_allele_names)
        
    def get_allowed_peptide_lengths_tests(self):
        ctlpanad = NetCTLpanAlleleData()
        
        # Filter only allele name
        allele_name = 'Patr-A0101'
        expected_lengths = [8, 9, 10, 11]
        lengths = ctlpanad.get_allowed_peptide_lengths(allele_name=allele_name)
        self.assertEqual(lengths, expected_lengths)
        
    def get_species_for_allele_name_tests(self):
        ctlpanad = NetCTLpanAlleleData()
        
        allele_name = 'HLA-A*02:01'
        species = ctlpanad.get_species_for_allele_name(allele_name)
        self.assertEqual(species, 'human')

        allele_name = 'Patr-A*0101'
        species = ctlpanad.get_species_for_allele_name(allele_name)
        self.assertEqual(species, 'chimpanzee')

class AlleleInfoFunctionTests(TestCase):
    def is_user_defined_allele_tests(self):
        test_cases = {
            'MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWI': True,
            'PEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQSMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRS': True,
            'WTAADTAAQITQRKWEAAREAEQRRAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGF': True,
            'HLA-A*11:01': False,
            'SLA-1-HB04': False,
            'RT1A': False,
        }

        for allele, expected in test_cases.items():
            result = is_user_defined_allele(allele)
            self.assertEqual(result, expected)
