from unittest import TestCase
import os
from .funcs import dcs
from . import KEEPSAKES
from itertools import permutations
class TestOutputFiles(TestCase):

    def setUp(self) -> None:
        [os.remove(f) for f in os.listdir() if f not in KEEPSAKES]
        return super().setUp()

    def tearDown(self) -> None:
        [os.remove(f) for f in os.listdir() if f not in KEEPSAKES]

    def dcs(self):
        """
        Wrapper to call dcs on constant path.
        """
        path = os.path.abspath("examples/NIMROD/NIMRODsilica/NIMRODmesoporoussilica.txt")
        return dcs(path)

    def make_combination(iterA, iterB):
        return [list(zip(perm, iterB) for perm in permutations(iterA, len(iterB)))]
    def testGudrunMakesGeneralOutputFiles(self):
        self.dcs()
        self.assertTrue("spike.dat" in os.listdir())
        self.assertTrue("deadtime.cor" in os.listdir())
        self.assertTrue("gudrun_run_par.dat" in os.listdir())
        self.assertTrue("gudrun_grp.dat" in os.listdir())
        self.assertTrue("gudrun_calib.dat" in os.listdir())
        self.assertTrue("gudrun_van_tcb.dat" in os.listdir())
        self.assertTrue("gudrun_sam_tcb.dat" in os.listdir())
        self.assertTrue("vanadium.soq" in os.listdir())
        self.assertTrue("vansmo.par" in os.listdir())
        self.assertTrue("ran1.dat" in os.listdir())
        self.assertTrue("polyfitcoeff.text" in os.listdir())
        self.assertTrue("runfactor_list.dat" in os.listdir())
        self.assertTrue("vanadiun.module" in os.listdir())
        self.assertTrue("vanadium_back.module" in os.listdir())
        self.assertTrue("sample_back.module" in os.listdir())
        self.assertTrue("sample.module" in os.listdir())

    def testGudrunMakesContainerSpecificOutputFiles(self):
        EXPECTED_PREFIXES = ["NIMROD00000772", "NIMROD00000772", "NIMROD00000772"]
        EXPECTED_SUFFIXES = [".mul01", ".mut01", ".rawmon", ".rawtrans", ".smomon", ".trans01", ".samrat" , ".transnofit01"]

        for file in self.make_combination(EXPECTED_PREFIXES, EXPECTED_SUFFIXES):
            self.assertTrue(file in os.listdir())

    def testGudrunMakesSampleSpecificOutputFiles(self):
        EXPECTED_PREFIXES = ["NIMROD00000769", "NIMROD00000770", "NIMROD00000771"]
        EXPECTED_SUFFIXES = [".abs01", ".bak", ".gr1", ".gr2", ".gud", ".mul01", ".mut01", ".pla01", ".rawmon", ".rawtrans", ".smomon", ".sub01", ".dcs01", ".mdcs01", ".mgor01", ".mdor01", ".msub01",
        ".msub01", ".chksum", ".samrat", ".transnofit01"]
        
        for file in self.make_combination(EXPECTED_PREFIXES, EXPECTED_SUFFIXES):
            self.assertTrue(file in os.listdir())

    def testGudrunMakesNormalisationSpecificOutputFiles(self):

        EXPECTED_PREFIXES = ["NIMROD00000767"]
        EXPECTED_SUFFIXES = [".mul01", ".mut01", ".pla01", ".rawmon", ".rawtrans", ".smomon", ".trans01", ".transnofit01", ".vanrat"]

        for file in self.make_combination(EXPECTED_PREFIXES, EXPECTED_SUFFIXES):
            self.assertTrue(file in os.listdir())
    
    def testGudrunMakesSampleBackgroundSpecificOutputFiles(self):

        EXPECTED_PREFIXES = ["NIMROD00000768"]
        EXPECTED_SUFFIXES = [".rawmon", ".rawtrans", ".smomon"]

        for file in self.make_combination(EXPECTED_PREFIXES, EXPECTED_SUFFIXES):
            self.assertTrue(file in os.listdir())

    # mint01
