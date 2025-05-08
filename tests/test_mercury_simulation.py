"""
Unit tests for the Mercury E-QFT simulation.
"""

import unittest
import numpy as np
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from mercury_eqft_simulation import (
    acceleration_newton,
    acceleration_gr_correction,
    acceleration_eqft_correction,
    G, c, M_SUN, r_s
)

class TestMercurySimulation(unittest.TestCase):
    """Test cases for the Mercury E-QFT simulation."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a simple test state at 1 AU from the sun
        self.AU = 149.6e9  # 1 AU in meters
        self.state = [self.AU, 0, 0, 0, 30000, 0]  # Position and velocity vector
        self.time = 0
    
    def test_acceleration_newton(self):
        """Test the Newtonian acceleration calculation."""
        accel = acceleration_newton(self.time, self.state)
        # Check direction (should point toward the sun)
        self.assertLess(accel[0], 0)
        self.assertAlmostEqual(accel[1], 0)
        self.assertAlmostEqual(accel[2], 0)
        
        # Check magnitude (F = GMm/r² -> a = GM/r²)
        expected_magnitude = G * M_SUN / (self.AU ** 2)
        actual_magnitude = np.sqrt(accel[0]**2 + accel[1]**2 + accel[2]**2)
        self.assertAlmostEqual(actual_magnitude, expected_magnitude, places=10)
    
    def test_acceleration_gr_correction(self):
        """Test the GR correction to acceleration."""
        accel_gr = acceleration_gr_correction(self.time, self.state)
        
        # GR correction should be much smaller than Newtonian acceleration
        accel_newton = acceleration_newton(self.time, self.state)
        newton_magnitude = np.sqrt(np.sum(np.array(accel_newton)**2))
        gr_magnitude = np.sqrt(np.sum(np.array(accel_gr)**2))
        
        # GR correction should be on the order of (r_s/r) times Newtonian
        ratio = gr_magnitude / newton_magnitude
        expected_ratio_order = r_s / self.AU
        self.assertLess(ratio, expected_ratio_order * 10)
    
    def test_acceleration_eqft_correction(self):
        """Test the E-QFT correction to acceleration."""
        accel_eqft = acceleration_eqft_correction(self.time, self.state)
        
        # E-QFT correction should be non-zero
        eqft_magnitude = np.sqrt(np.sum(np.array(accel_eqft)**2))
        self.assertGreater(eqft_magnitude, 0)
        
        # E-QFT correction should be smaller than GR correction
        accel_gr = acceleration_gr_correction(self.time, self.state)
        gr_magnitude = np.sqrt(np.sum(np.array(accel_gr)**2))
        self.assertLess(eqft_magnitude, gr_magnitude)

if __name__ == '__main__':
    unittest.main()