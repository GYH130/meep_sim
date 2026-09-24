"""Offline deadline/cap invariants; never imports Meep."""
import unittest
from run_k import cap_seconds, SOLVER_DEADLINE, RENTAL_EXPIRY

class PolicyTests(unittest.TestCase):
    def test_caps(self):
        self.assertEqual(cap_seconds({'single_cap_seconds':14400}, 30000), 14400)
        self.assertEqual(cap_seconds({'single_cap_seconds':21600}, 5000), 5000)
        self.assertEqual(cap_seconds({'single_cap_seconds':99999}, 90000), 21600)
        self.assertEqual(cap_seconds({'single_cap_seconds':21600}, -1), 0)

    def test_backup_reserve(self):
        self.assertEqual(RENTAL_EXPIRY-SOLVER_DEADLINE, 83*60)

if __name__ == '__main__':
    unittest.main()
