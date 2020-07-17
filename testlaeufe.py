#!/Users/alexfinck/anaconda3/bin/python
# -*- coding: utf-8 -*-

import show_alignment
import unittest
import random


class ShowAlignment_Test(unittest.TestCase):
    def setUp(self):
        self.aligner = show_alignment.ShowAlignment()
        
    def testCase1(self):
        print("\ntestCase1> simple example")
        # Zeile 1 vom Zettel - dient nur zu Testzwecken
        obereZeile = "GGTCTAAG_TA"
        # zeile 2 dient nur zu Testzwecken
        untereZeile = "GGACTAAGGTA"
        self.assertTrue(
            self.aligner.showAlignment(obereZeile, untereZeile)
        )
       
    def testCase2(self):
        print("\n\ntestCase2> only gaps")
        obereZeile = "___________"
        untereZeile = "GGACTAAGGTA"
        self.assertTrue(
            self.aligner.showAlignment(obereZeile, untereZeile)
        )
        
    def testCase3(self):
        print("\n\ntestCase3> input error with aa in upper line given a dna alphabet")
        obereZeile = "GGTBTAAG_TA"
        untereZeile = "GGACTAAGGTA"
        with self.assertRaises(ValueError):
            self.aligner.showAlignment(obereZeile, untereZeile)

    def testCase4(self):
        print("\n\ntestCase4> input error with aa in lower line given a dna alphabet")
        obereZeile = "GGACTAAGGTA"
        untereZeile = "GGTBTAAG_TA"
        with self.assertRaises(ValueError):
            self.aligner.showAlignment(obereZeile, untereZeile)

    def testCase5(self):
        print("\n\ntestCase5> upper and lower sequence have different input length")
        obereZeile =  "__GGACTAAGGTA"
        untereZeile = "GGTBTAAG_TA"
        with self.assertRaises(ValueError):
            self.aligner.showAlignment(obereZeile, untereZeile)

    def testCase6(self):
        print("\n\ntestCase6> upper sequence is made of blanks")
        obereZeile = "           "
        untereZeile = "XGACTAAGGTA"
        with self.assertRaises(ValueError):
            self.aligner.showAlignment(obereZeile, untereZeile)

    def testCase7(self):
        print("\n\ntestCase7> invalid input in upper sequence with protein alphabet")
        obereZeile  = "JGMBTAAG_TA"
        untereZeile = "XGACTAAGGTA"
        self.aligner.aligntIsDna = False
        with self.assertRaises(ValueError):
            self.aligner.showAlignment(obereZeile, untereZeile)

    def testCase8(self):
        print("\n\ntestCase8> multiline input with dna alphabet")
        obereZeile  = "GGACTAAG_TAGGACTAAG_TAGGACTAAG_TAGGACTAAG_TAGGAC" + \
            "_TAGGACTAAG_TAGGACTAAG_TAGGA"
        untereZeile = "GACTAAGGTAGACTAAGGTAGACTAAGGTAGACTAAGGTAGACTAAGG" + \
            "CTAAGGTAGACTAAGGTAGACTAAGGTA"
        self.assertTrue(
            self.aligner.showAlignment(obereZeile, untereZeile)
        )

    def testCase9(self):
        print("\n\ntestCase9> multiline input with dna alphabet")
        self.aligner.aligntIsDna = False
        obereZeile  = "GGACTAAG_TAGGACTAAG_TAGGACTAAG_TAGGACTAAG_TAGGACTA" + \
            "_TAGGACTAAG_TAGGACTAAG_TAGGA"
        untereZeile = "GACTAAGGTAGACTAAGGTAGACTAAGGTAGACTAAGGTAGACTAAGGTA" + \
            "CTAAGGTAGACTAAGGTAGACTAAGGTA"
        self.assertTrue(
            self.aligner.showAlignment(obereZeile, untereZeile)
        )

    def testCase10(self):
        print("\n\ntestCase10> long input with protein alphabet")
        self.aligner.aligntIsDna = False
        obereZeile  = ''.join(random.choices(
            list(self.aligner.PROTEIN_ALPHABET), k=600
        ))
        untereZeile = ''.join(random.choices(
            list(self.aligner.PROTEIN_ALPHABET), k=600
        ))
        self.assertTrue(
            self.aligner.showAlignment(obereZeile, untereZeile)
        )
            
    def testCase11(self):
        print("\n\ntestCase11> long input with dna alphabet")
        self.aligner.aligntIsDna = False
        obereZeile  = ''.join(random.choices(
            list(self.aligner.DNA_ALPHABET), k=600
        ))
        untereZeile = ''.join(random.choices(
            list(self.aligner.DNA_ALPHABET), k=600
        ))
        self.assertTrue(
            self.aligner.showAlignment(obereZeile, untereZeile)
        )
        
    def testCase12(self):
        print("\n\ntestCase12> lowercase input")
        self.aligner.aligntIsDna = False
        obereZeile  = ''.join(random.choices(
            list(str.lower(self.aligner.DNA_ALPHABET)), k=60
        ))
        untereZeile = ''.join(random.choices(
            list(str.lower(self.aligner.DNA_ALPHABET)), k=60
        ))
        self.assertTrue(
            self.aligner.showAlignment(obereZeile, untereZeile)
        )

        
def suite():
    suite = unittest.TestSuite()
    suite.addTest(ShowAlignment_Test('testCase1'))
    suite.addTest(ShowAlignment_Test('testCase2'))
    suite.addTest(ShowAlignment_Test('testCase3'))
    suite.addTest(ShowAlignment_Test('testCase4'))
    suite.addTest(ShowAlignment_Test('testCase5'))
    suite.addTest(ShowAlignment_Test('testCase6'))
    suite.addTest(ShowAlignment_Test('testCase7'))
    suite.addTest(ShowAlignment_Test('testCase8'))
    suite.addTest(ShowAlignment_Test('testCase9'))
    suite.addTest(ShowAlignment_Test('testCase10'))
    suite.addTest(ShowAlignment_Test('testCase11'))
    suite.addTest(ShowAlignment_Test('testCase12'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
