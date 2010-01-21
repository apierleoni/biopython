#!/usr/bin/env python
"""Test for the Uniprot parser on Uniprot XML files.
"""
import os
import copy
import unittest

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from seq_tests_common import compare_records

class TestUniprot(unittest.TestCase):

    def test_uni001(self):
        "Parsing Uniprot file uni001"
        filename = 'uni001'
        # test the record parser

        datafile = os.path.join('Uniprot', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "uniprot")
        test_handle.close()

        self.assert_(isinstance(seq_record, SeqRecord))

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(seq_record.id, "Q91G55")
        self.assertEqual(seq_record.name, "043L_IIV6")
        self.assertEqual(seq_record.description, "Uncharacterized protein 043L")
        self.assertEqual(repr(seq_record.seq), "Seq('MDLINNKLNIEIQKFCLDLEKKYNINYNNLIDLWFNKESTERLIKCEVNLENKI...IPI', ProteinAlphabet())")

        # self.assertEqual(seq_record.accessions, ['Q91G55']) #seq_record.accessions does not exist
        # self.assertEqual(seq_record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        # self.assertEqual(record.seqinfo, (348, 39676, '75818910'))
    
        self.assertEqual(len(seq_record.features), 1)
        self.assertEqual(repr(seq_record.features[0]), "SeqFeature(FeatureLocation(ExactPosition(0),ExactPosition(115)), type='chain', id='PRO_0000377969')")

        self.assertEqual(len(seq_record.annotations['references']), 2)
        self.assertEqual(seq_record.annotations['references'][0].authors, 'Jakob N.J., Mueller K., Bahr U., Darai G.')
        self.assertEqual(seq_record.annotations['references'][0].title, 'Analysis of the first complete DNA sequence of an invertebrate iridovirus: coding strategy of the genome of Chilo iridescent virus.')
        self.assertEqual(seq_record.annotations['references'][0].journal, 'Virology 286:182-196(2001)')
        self.assertEqual(seq_record.annotations['references'][0].comment, 'journal article | 2001 | Scope: NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA] | ')

        self.assertEqual(len(seq_record.dbxrefs), 11)
        self.assertEqual(seq_record.dbxrefs[0], 'DOI:10.1006/viro.2001.0963')

        self.assertEqual(seq_record.annotations['sequence_length'], 116)
        self.assertEqual(seq_record.annotations['sequence_checksum'], '4A29B35FB716523C')
        self.assertEqual(seq_record.annotations['modified'], '2009-07-07')
        self.assertEqual(seq_record.annotations['accessions'], ['Q91G55'])
        self.assertEqual(seq_record.annotations['taxonomy'], ['Viruses', 'dsDNA viruses, no RNA stage', 'Iridoviridae', 'Iridovirus'])
        self.assertEqual(seq_record.annotations['sequence_mass'], 13673)
        self.assertEqual(seq_record.annotations['dataset'], 'Swiss-Prot')
        self.assertEqual(seq_record.annotations['gene_name_ORF'], ['IIV6-043L'])
        self.assertEqual(seq_record.annotations['version'], 21)
        self.assertEqual(seq_record.annotations['sequence_modified'], '2001-12-01')
        self.assertEqual(seq_record.annotations['keywords'], ['Complete proteome', 'Virus reference strain'])
        self.assertEqual(seq_record.annotations['organismHost_name'], ['Acheta domesticus', 'House cricket', 'Chilo suppressalis', 'striped riceborer', 'Gryllus bimaculatus', 'Two-spotted cricket', 'Gryllus campestris', 'Spodoptera frugiperda', 'Fall armyworm'])
        self.assertEqual(seq_record.annotations['created'], '2009-06-16')
        self.assertEqual(seq_record.annotations['organism_name'], ['Chilo iridescent virus'])
        self.assertEqual(seq_record.annotations['organism'], 'Invertebrate iridescent virus 6 (IIV-6)')
        self.assertEqual(seq_record.annotations['recommendedName_fullName'], ['Uncharacterized protein 043L'])
        self.assertEqual(seq_record.annotations['sequence_version'], 1)
        self.assertEqual(seq_record.annotations['proteinExistence'], ['Predicted'])

    def test_Q13639(self):
        up_filename = "Uniprot/Q13639.xml"
        sp_filename = "SwissProt/Q13639.txt"
        up_records = list(SeqIO.parse(open(up_filename),'uniprot'))
        sp_records = list(SeqIO.parse(open(sp_filename),'swiss'))
        compare_records(sp_records, up_records)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
