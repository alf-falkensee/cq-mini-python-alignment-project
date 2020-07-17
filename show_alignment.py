#!/Users/alexfinck/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""This script extends an alignment consisting of upper-seq and
lower-seq with a quality line shown between upper-seq and lower-seq on
a terminal or console having at least 80 columns and configured to use
a fixed/monospaced font.

The script takes as main optional argument whether or not an input is
dna or proteins; this is an authoritative argument as the script does
not try to guess in place of the user; default is dna. The script
returns a message to the console consisting in three sections. For
each run the user is informed about pre-requisites to a succesful
usage in a first section. Then the input is parsed for consistency and
a corresponding section is printed to the console. The last section
contains the input alignement accompanied by a quality line in the
form of 3 lines repeatedly displayed and cut at default 60 columns
until the length of the alignment is reached.

"""
# "Skelett" des Programs zum Zeigen eines Alignments mit einer
# Qualitaetszeile Haltet euch bitte an diese Struktur - dort wo das
# "pass" steht muss euer Quelltext kommen. Das pass muss dazu
# geloescht werden.

# Autor: Alex Finck
# Datum der letzten Aenderung: 09.07.2020
#
# usage examples from the command line:
# > python show_alignment.py "AACTG_GTCAT" "AGTCAA_CTGA"
# > python show_alignment.py -iprotein "ACTG_GTCA" "GTCAA_CTG"


from Bio.SubsMat.MatrixInfo import pam30
import argparse


class ShowAlignment:
    def __init__(self, aligntIsDna=True):
        # Liste die angibt wie gut ein Match ist. Siehe Aufgabenzettel -
        self.OUTPUT_WIDTH_IN_COLS = 60
        self.qualitaetsListe = ["AA", "GG", "CC", "TT", "CT", "TC", "AG", "GA",
                                "CA", "AC", "CG", "GC", "TA", "AT", "TG", "GT"]
        self.qualitaetsListeProteins = self.__qualitaetsListeProteins()
        # beware: the following variable has coupling with the above
        # _qualitaetsListeProteins function
        self.PROTEIN_ALPHABET = 'ABCDEFGHIKLMNPQRSTVWXYZ_'
        self.DNA_ALPHABET = "ACTG_"
        self.INPUT_GAP_ZEICHEN = "_"
        self.DNA_EXAKTER_MATCH = range(0, 4)
        self.DNA_GUTER_MATCH = range(4, 8)
        self.DNA_KEIN_GUTER_MATCH = range(8, 16)
        # beware: the following 3 variables have coupling with the above
        # _qualitaetsListeProteins function
        self.AA_EXAKTER_MATCH = range(0, 19)
        self.AA_GUTER_MATCH = range(19, 60)
        self.AA_KEIN_GUTER_MATCH = range(60, 529)
        
        self.EXAKTER_MATCH_ZEICHEN = "|"
        self.GUTER_MATCH_ZEICHEN = ":"
        self.KEIN_GUTER_MATCH_ZEICHEN = "."
        
        self.QUAL_GAP_ZEICHEN = " "
        
        # flags, the values of which are given back to the user 
        self.VALID_DNA_OR_PROTEIN = "valid dna|protein"
        self.INVALID_PROTEIN = "invalid protein"
        self.INVALID_DNA = "invalid dna"
        # deductive and authoritative flag that gives to the script the prior
        # information about 
        self.aligntIsDna = True
    
    # string[] __qualitaetsListeProteins()
    def __qualitaetsListeProteins(self):
        """Private function building and returning a quality list analog to
        the quality list for dna nucleotides, but based upon the PAM30
        Matrix; associated quaity ranges are defined in AA_EXAKTER_MATCH,
        AA_GUTER_MATCH, AA_KEIN_GUTER_MATCH and correspond for
        AA_GUTER_MATCH to remaining positve scores after removal of exact
        matches and for AA_KEIN_GUTER_MATCH to negative scores,
        respectively        
        """
        rv = []
        pam30_sortierbar = {}
        for key in pam30.keys():
            pam30_sortierbar[str(pam30[key]) + ";" + ''.join(key)] = pam30[key]
            if key[0] != key[1]:
                pam30_sortierbar[
                    str(pam30[key]) + ";" + ''.join((key[1], key[0]))
                ] = pam30[key]
        sorted_keys = list(pam30_sortierbar.keys())
        sorted_keys.sort(key=lambda k: int(k.split(";")[0]), reverse=True)
        # debugging kept for historical reasons
        #    for key in iter(sorted_keys):
        #        print(key.split(";")[1] + " has score " + str(pam30_sortierbar[key]))
        for key in iter(sorted_keys):
            rv.append(key.split(";")[1])
        return(rv)

    # string getQuality(string obereZeile, string untereZeile)
    def getQuality(self, obereZeile, untereZeile):
        """Function that returns in the form of a string a quality of an
        alignment consisting of two input sequences of dna or proteins.
        The quality depnds on the prior input of the user given by the
        aligntIsDna Flag. Quality for dna sequence pairs further depends
        upon the list 'qualitaetsListe' and for amino acid sequences upon
        the list 'qualitaetsListeProteins'.
        """
        qualitaetsZeile = ""
        if self.aligntIsDna:
            _exakter_match_list = self.DNA_EXAKTER_MATCH
            _guter_match_list = self.DNA_GUTER_MATCH
            _kein_guter_match_list = self.DNA_KEIN_GUTER_MATCH
            _qualitaetsListe = self.qualitaetsListe
        else:
            _exakter_match_list = self.AA_EXAKTER_MATCH
            _guter_match_list = self.AA_GUTER_MATCH
            _kein_guter_match_list = self.AA_KEIN_GUTER_MATCH
            _qualitaetsListe = self.qualitaetsListeProteins
            
        for i in range(len(obereZeile)):
            if (
                obereZeile[i] == self.INPUT_GAP_ZEICHEN or
                untereZeile[i] == self.INPUT_GAP_ZEICHEN 
            ):
                qualitaetsZeile += self.QUAL_GAP_ZEICHEN
            else:
                currentResiduePair = str.upper(obereZeile[i] + untereZeile[i])
                #            print(currentResiduePair)
                indexOfPair = _qualitaetsListe.index(currentResiduePair)
                if indexOfPair in _exakter_match_list:
                    qualitaetsZeile += self.EXAKTER_MATCH_ZEICHEN
                if indexOfPair in _guter_match_list:
                    qualitaetsZeile += self.GUTER_MATCH_ZEICHEN
                if indexOfPair in _kein_guter_match_list:
                    qualitaetsZeile += self.KEIN_GUTER_MATCH_ZEICHEN
        return(qualitaetsZeile)
   
    # bool|ValueError showAlignment(cls, string zeile1, string zeile2)
    def showAlignment(self, zeile1, zeile2):
        """Function that processes an existing alignment of dna or proteins
        into a console output projection to a quality space determined by
        the function getQuality and typically consisting of zeile1 on the
        top, zeile2 at the bottom and a quality string in between. The
        console output is also separated in two sections. The first
        section is giving a feedback to the user about the consistency of
        the input alignment. In case of a consistent alignment and user
        choice (between dna and protein) a second section is displayed
        showing the alignment together with its quality within 60 columns
        of a properly (fixed fonts, more than 60 columns) configured
        console. In case of a succesful output to the console the function
        returns True. In case the consistency of the input is falsified,
        the first section is gracefully given back to the user, but
        processing of the input is interrupted by a ValueError exception.
    
        """

        if (self.inputCheckpoint(zeile1, zeile2)):
            # get the quality
            quality_zeile = self.getQuality(zeile1, zeile2)
            start_index = 0
            cutter_index = self.OUTPUT_WIDTH_IN_COLS
            while (start_index < len(quality_zeile)):
                print(zeile1[start_index:cutter_index])
                print(quality_zeile[start_index:cutter_index])
                print(zeile2[start_index:cutter_index])
                start_index = cutter_index
                targeted_end_index = cutter_index + self.OUTPUT_WIDTH_IN_COLS
                if targeted_end_index <= len(quality_zeile):
                    cutter_index = targeted_end_index
                else:
                    cutter_index = len(quality_zeile)
        return True

    # {residueIndex : int, residue : char, recognizedAlphabet : string} getValidityOfResiduesInSequence(string seq)
    def getValidityOfResiduesInSequence(self, seq):
        """Function returning the consistency of an individual input sequence
        as a dictionary containing in the inconsistent case the residue
        location and value of the first inconsistency and values
        confirming the validity of the input sequence otherwise.
        """
        seqList = list(seq)
        aSpotted_Index = -1
        aSpotted_residue = ""
        if self.aligntIsDna:
            _alphabet = self.DNA_ALPHABET
        else:
            _alphabet = self.PROTEIN_ALPHABET
        # iterate over the sequence given the prior knowldege of the user
        for i in range(len(seqList)):
            residue = seqList[i]
            if str.upper(residue) not in list(_alphabet):
                aSpotted_Index = i
                aSpotted_residue = residue
                break
        rv = {
            "residueIndex": aSpotted_Index,
            "residue": aSpotted_residue,
            "recognizedAlphabet": self.VALID_DNA_OR_PROTEIN
        }
        if (aSpotted_residue != ""):
            if self.aligntIsDna:
                rv["recognizedAlphabet"] = self.INVALID_DNA
            else:
                rv["recognizedAlphabet"] = self.INVALID_PROTEIN
        return(rv)

    # bool|ValueError inputCheckpoint(string obereZeile, string untereZeile)
    def inputCheckpoint(self, obereZeile, untereZeile):
        """Function checking the consistency of an alignment and generating
        output of the first section of showAlignment in its behalf. If an
        inconsistency is detected information about reasons for stopping
        further processing is given back to the user and a ValueError is
        raised. In case no inconsistency is found a summary report is also
        generated and the function returns True. The function accepts 1)
        only equal length for obereZeile, untereZeile 2) only the input
        alphabets + INPUT_GAP_ZEICHEN ("_")
        """
        rv = True
        # 1) only equal length for obereZeile, untereZeile
        if (len(obereZeile) != len(untereZeile)):
            print("============================================================")
            print("input sequences do not have the same length")
            print("============================================================")
            raise ValueError("Input sequences of different lengths")
        
        # 2) only the input alphabets + INPUT_GAP_ZEICHEN ("_")
        validityInObereZeile = self.getValidityOfResiduesInSequence(obereZeile)
        validityInUntereZeile = self.getValidityOfResiduesInSequence(untereZeile)
        if (
            validityInObereZeile["recognizedAlphabet"] == self.VALID_DNA_OR_PROTEIN and
            validityInUntereZeile["recognizedAlphabet"] == self.VALID_DNA_OR_PROTEIN
        ):
            print("============================================================")
            print("input is recognized as: " + self.VALID_DNA_OR_PROTEIN)
            _input_type = "dna"
            if not self.aligntIsDna:
                _input_type = "protein"
            print("input is now further processed as: " + _input_type)
            print("============================================================ ")
        else:
            print("============================================================ ")
            if (
                validityInObereZeile["recognizedAlphabet"] in
                [self.INVALID_DNA, self.INVALID_PROTEIN]
            ):
                print(
                    "upper sequence is recognized as: " +
                    validityInObereZeile["recognizedAlphabet"]
                )
                print(
                    "character number {} with value '{}' could not be parsed".
                    format(
                        validityInObereZeile["residueIndex"] + 1,
                        validityInObereZeile["residue"]
                    )
                )
            if (
                validityInUntereZeile["recognizedAlphabet"] in
                [self.INVALID_DNA, self.INVALID_PROTEIN]
            ):
                print(
                    "lower sequence is recognized as: " +
                    validityInUntereZeile["recognizedAlphabet"]
                )
                print(
                    "character number {} with value '{}' could not be parsed".
                    format(
                        validityInUntereZeile["residueIndex"] + 1,
                        validityInUntereZeile["residue"]
                    )
                )
            print("============================================================ ")
            raise ValueError("Input outside of chosen alphabet.")
        return(rv)

    # None informUserAboutPrerequisites()
    @staticmethod
    def informUserAboutPrerequisites():
        """A pure side-effect function informing the user about proper use of
        the script since no extra-care is taken by the script in order to
        enforce a proper console configuration
    
        """
        print("============================================================ ")
        print("make sure your terminal is set to use fixed/monospaced fonts ")
        print("and displays a minimum of 60 columns")
        print("============================================================ ")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        metavar='upper_seq', type=str, nargs=1, dest='upper_seq', 
        help='a first sequence of an input alignment displayed' + 
        ' on the top of the quality line'
    )
    parser.add_argument(
        metavar='lower_seq', type=str, nargs=1, dest='lower_seq', 
        help='a second sequence of an input ' +
        'alignment displayed below the quality line'
    )
    parser.add_argument(
        '-i', '--input_type', default="dna", type=str, 
        help='dna|protein defaults to dna'
    )
    args = parser.parse_args()
#
#    print(args.input_type)
#    print(args.upper_seq)
#    print(args.lower_seq)
    obereZeile = args.upper_seq[0]
    untereZeile = args.lower_seq[0]
    if (args.input_type == "dna"):
        aligntIsDna = True
    elif (args.input_type == "protein"):
        aligntIsDna = False
    ShowAlignment.informUserAboutPrerequisites()
    a = ShowAlignment(aligntIsDna)
    a.showAlignment(obereZeile, untereZeile)
