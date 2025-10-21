#!/usr/bin/env python3
"""
NiTi Nanoparticle Generator with Defects
Author: PrudhviNallagatla
Date: 2025-10-21

Generates spherical nitinol (NiTi) nanoparticles with various defect types
Outputs XYZ file for LAMMPS simulations
"""

import numpy as np
import argparse
from dataclasses import dataclass
from typing import List, Tuple, Optional
import sys


@dataclass
class Atom:
    """Class to represent an atom"""

    element: str
    position: np.ndarray
    atom_id: int


class NiTiNanoparticle:
    """
    Generate spherical NiTi nanoparticles with B2 (austenite) structure
    """

    def __init__(
        self,
        diameter_nm: float,
        ni_percent: float = 50.0,
        lattice_param: float = 3.015,
        shape: str = "sphere",
        aspect_ratio: float = 1.0,
    ):
        """
        Initialize nanoparticle generator

        Args:
            diameter_nm: Desired diameter in nanometers
            ni_percent: Percentage of Ni atoms (0-100, default: 50.0 for equiatomic)
            lattice_param: B2 lattice parameter in Angstroms (default: 3.015 Å)
            shape: Particle shape - "sphere", "ellipsoid", "faceted", or "rough"
            aspect_ratio: For ellipsoid shape (default: 1.0)
        """
        self.diameter_nm = diameter_nm
        self.diameter_angstrom = diameter_nm * 10.0
        self.radius = self.diameter_angstrom / 2.0
        self.lattice_param = lattice_param
        self.ni_percent = ni_percent
        self.ti_percent = 100.0 - ni_percent
        self.atoms = []
        self.shape = shape
        self.aspect_ratio = aspect_ratio

        # Set up shape parameters
        if shape == "ellipsoid":
            self.ellipsoid_axes = np.array(
                [self.radius, self.radius * aspect_ratio, self.radius / aspect_ratio]
            )

        # Calculate supercell size needed
        self.n_cells = int(np.ceil(self.diameter_angstrom / self.lattice_param)) + 2

        print(f"Generating {diameter_nm} nm NiTi nanoparticle")
        print(f"Shape: {shape.capitalize()}")
        if shape == "ellipsoid":
            print(f"Aspect ratio: {aspect_ratio:.2f}")
        print(f"Target composition: Ni{ni_percent:.1f}Ti{self.ti_percent:.1f}")
        print(f"Lattice parameter: {lattice_param} Å")
        print(f"Supercell size: {self.n_cells}x{self.n_cells}x{self.n_cells}")

    def set_shape_type(self, shape: str = "sphere", aspect_ratio: float = 1.0):
        """
        Set the nanoparticle shape

        Args:
            shape: "sphere", "ellipsoid", "faceted", or "rough"
            aspect_ratio: For ellipsoid, ratio of axes (default: 1.0)
        """
        self.shape = shape
        self.aspect_ratio = aspect_ratio

        if shape == "ellipsoid":
            # Random orientation for ellipsoid
            self.ellipsoid_axes = np.array(
                [self.radius, self.radius * aspect_ratio, self.radius / aspect_ratio]
            )
            print(f"Shape: Ellipsoid with aspect ratio {aspect_ratio:.2f}")
        elif shape == "faceted":
            print("Shape: Faceted (Wulff-like)")
        elif shape == "rough":
            print("Shape: Rough surface")
        else:
            print("Shape: Spherical")

    def _is_inside_shape(self, position: np.ndarray) -> bool:
        """
        Check if position is inside the nanoparticle shape

        Args:
            position: 3D position vector

        Returns:
            True if inside, False otherwise
        """
        if self.shape == "sphere":
            return np.linalg.norm(position) <= self.radius

        elif self.shape == "ellipsoid":
            # Ellipsoid equation: (x/a)^2 + (y/b)^2 + (z/c)^2 <= 1
            normalized = position / self.ellipsoid_axes
            return np.dot(normalized, normalized) <= 1.0

        elif self.shape == "faceted":
            # Approximate Wulff construction for B2 structure
            # Favor {100} and {110} planes
            distance = np.linalg.norm(position)
            if distance > self.radius * 1.1:
                return False

            # Apply faceting: truncate at crystal planes
            abs_pos = np.abs(position)

            # {100} planes - truncate at faces
            max_coord = np.max(abs_pos)
            if max_coord > self.radius * 0.85:
                return distance <= self.radius * 0.95

            # {110} planes - truncate at edges
            sum_two_largest = np.sum(np.sort(abs_pos)[-2:])
            if sum_two_largest > self.radius * 1.2:
                return False

            return distance <= self.radius

        elif self.shape == "rough":
            # Add MILD surface roughness (not extreme!)
            distance = np.linalg.norm(position)

            # First check if clearly outside
            if distance > self.radius * 1.05:
                return False

            # Calculate roughness based on spherical harmonics
            theta = np.arctan2(position[1], position[0])
            phi = np.arccos(np.clip(position[2] / (distance + 1e-10), -1, 1))

            # REDUCED roughness amplitudes (was 0.03, 0.02, 0.015)
            roughness = (
                0.01 * np.sin(3 * theta) * np.cos(2 * phi)
                + 0.008 * np.cos(5 * theta) * np.sin(3 * phi)
                + 0.006 * np.sin(4 * phi)
            )

            effective_radius = self.radius * (1.0 + roughness)
            return distance <= effective_radius

        elif self.shape == "blob":
            # Irregular blob shape using spherical harmonics
            # This mimics rapid solidification from EDM liquid droplet
            distance = np.linalg.norm(position)

            # Calculate spherical coordinates
            if distance < 1e-10:
                return True  # Center is always inside

            theta = np.arctan2(position[1], position[0])
            phi = np.arccos(np.clip(position[2] / distance, -1, 1))

            # Use low-order spherical harmonics for smooth, organic shape
            # Coefficients chosen to create ~15-25% radius variation
            radius_factor = 1.0
            radius_factor += 0.15 * np.sin(2 * theta) * np.cos(phi)
            radius_factor += 0.12 * np.cos(3 * theta) * np.sin(2 * phi)
            radius_factor += 0.10 * np.sin(theta) * np.cos(3 * phi)
            radius_factor += 0.08 * np.cos(2 * theta) * np.sin(phi)
            radius_factor += 0.06 * np.sin(4 * phi)

            # Ensure radius_factor stays positive (convex shape)
            radius_factor = max(0.75, min(1.25, radius_factor))

            effective_radius = self.radius * radius_factor
            return distance <= effective_radius

        return np.linalg.norm(position) <= self.radius

    def generate_b2_structure(self):
        """
        Generate B2 (CsCl-type) crystal structure for NiTi austenite
        B2 structure: Simple cubic with Ni at corners, Ti at body center (or vice versa)
        Then adjust composition to match desired Ni:Ti ratio
        """
        print("Generating B2 austenite structure...")

        atom_id = 1
        a = self.lattice_param

        # Calculate center of the supercell
        center = np.array([self.n_cells * a / 2.0] * 3)

        # First, generate perfect B2 structure
        temp_atoms = []

        for i in range(self.n_cells):
            for j in range(self.n_cells):
                for k in range(self.n_cells):
                    # Base position
                    base_pos = np.array([i * a, j * a, k * a])

                    # Ni at corner (0,0,0)
                    ni_pos = base_pos - center  # Relative to center
                    if self._is_inside_shape(ni_pos):
                        temp_atoms.append(
                            {"element": "Ni", "position": (ni_pos + center).copy()}
                        )

                    # Ti at body center (0.5, 0.5, 0.5)
                    ti_pos = base_pos + np.array([a / 2, a / 2, a / 2]) - center
                    if self._is_inside_shape(ti_pos):
                        temp_atoms.append(
                            {"element": "Ti", "position": (ti_pos + center).copy()}
                        )

        print(f"Generated {len(temp_atoms)} atoms in perfect B2 spherical region")

        # Now adjust composition to match target Ni percentage
        self._adjust_composition(temp_atoms)

        # Create final atom list with IDs
        for i, atom_data in enumerate(temp_atoms, 1):
            self.atoms.append(Atom(atom_data["element"], atom_data["position"], i))

        # Recenter particle at origin
        self._recenter_particle()

        # Check composition
        self._check_composition()

    def _adjust_composition(self, temp_atoms):
        """
        Adjust composition to match target Ni percentage
        by converting some atoms from one type to another
        """
        # Count current atoms
        current_ni = sum(1 for a in temp_atoms if a["element"] == "Ni")
        current_ti = sum(1 for a in temp_atoms if a["element"] == "Ti")
        total = len(temp_atoms)

        target_ni = int(total * self.ni_percent / 100.0)
        target_ti = total - target_ni

        diff_ni = target_ni - current_ni

        print(f"\nAdjusting composition:")
        print(f"  Current: Ni={current_ni}, Ti={current_ti}")
        print(f"  Target:  Ni={target_ni}, Ti={target_ti}")
        print(f"  Need to convert {abs(diff_ni)} atoms")

        if diff_ni > 0:
            # Need more Ni, convert Ti -> Ni
            ti_indices = [i for i, a in enumerate(temp_atoms) if a["element"] == "Ti"]
            if len(ti_indices) < diff_ni:
                print(f"  WARNING: Not enough Ti atoms to convert!")
                diff_ni = len(ti_indices)

            indices_to_convert = np.random.choice(ti_indices, diff_ni, replace=False)
            for idx in indices_to_convert:
                temp_atoms[idx]["element"] = "Ni"

        elif diff_ni < 0:
            # Need more Ti, convert Ni -> Ti
            ni_indices = [i for i, a in enumerate(temp_atoms) if a["element"] == "Ni"]
            if len(ni_indices) < abs(diff_ni):
                print(f"  WARNING: Not enough Ni atoms to convert!")
                diff_ni = -len(ni_indices)

            indices_to_convert = np.random.choice(
                ni_indices, abs(diff_ni), replace=False
            )
            for idx in indices_to_convert:
                temp_atoms[idx]["element"] = "Ti"

    def _recenter_particle(self):
        """Move center of mass to origin"""
        if not self.atoms:
            return

        com = np.mean([atom.position for atom in self.atoms], axis=0)
        for atom in self.atoms:
            atom.position -= com

        print("Particle recentered at origin")

    def _check_composition(self):
        """Check and report Ni:Ti ratio"""
        n_ni = sum(1 for atom in self.atoms if atom.element == "Ni")
        n_ti = sum(1 for atom in self.atoms if atom.element == "Ti")
        total = n_ni + n_ti

        actual_ni_percent = 100.0 * n_ni / total if total > 0 else 0
        actual_ti_percent = 100.0 * n_ti / total if total > 0 else 0

        print(f"\nFinal Composition:")
        print(f"  Ni: {n_ni} atoms ({actual_ni_percent:.2f}%)")
        print(f"  Ti: {n_ti} atoms ({actual_ti_percent:.2f}%)")
        print(f"  Total atoms: {total}")
        print(f"  Composition: Ni{actual_ni_percent:.1f}Ti{actual_ti_percent:.1f}")

    def add_vacancies(self, concentration: float = 0.01):
        """
        Add vacancy defects by removing random atoms

        Args:
            concentration: Fraction of atoms to remove (default: 0.01 = 1%)
        """
        if concentration <= 0:
            return

        n_initial = len(self.atoms)
        n_remove = int(n_initial * concentration)

        print(
            f"\nAdding vacancies: removing {n_remove} atoms ({100*concentration:.1f}%)"
        )

        # Randomly select atoms to remove
        indices_to_remove = np.random.choice(len(self.atoms), n_remove, replace=False)
        indices_to_remove = sorted(indices_to_remove, reverse=True)

        for idx in indices_to_remove:
            del self.atoms[idx]

        # Reassign atom IDs
        for i, atom in enumerate(self.atoms, 1):
            atom.atom_id = i

        print(f"Vacancies created: {n_initial} -> {len(self.atoms)} atoms")

    def add_antisite_defects(self, concentration: float = 0.02):
        """
        Add antisite defects by swapping Ni and Ti atoms

        Args:
            concentration: Fraction of atoms to swap (default: 0.02 = 2%)
        """
        if concentration <= 0:
            return

        print(f"\nAdding antisite defects: {100*concentration:.1f}% swaps")

        # Get indices of Ni and Ti atoms
        ni_indices = [i for i, atom in enumerate(self.atoms) if atom.element == "Ni"]
        ti_indices = [i for i, atom in enumerate(self.atoms) if atom.element == "Ti"]

        # Number of pairs to swap
        n_swaps = int(len(self.atoms) * concentration / 2)
        n_swaps = min(n_swaps, len(ni_indices), len(ti_indices))

        # Randomly select pairs to swap
        ni_to_swap = np.random.choice(ni_indices, n_swaps, replace=False)
        ti_to_swap = np.random.choice(ti_indices, n_swaps, replace=False)

        # Perform swaps
        for ni_idx, ti_idx in zip(ni_to_swap, ti_to_swap):
            self.atoms[ni_idx].element = "Ti"
            self.atoms[ti_idx].element = "Ni"

        print(f"Created {n_swaps} antisite pairs ({2*n_swaps} atoms affected)")

    def add_surface_vacancies(self, concentration: float = 0.05):
        """
        Add vacancies preferentially at surface

        Args:
            concentration: Fraction of surface atoms to remove (default: 0.05 = 5%)
        """
        if concentration <= 0:
            return

        print(f"\nAdding surface vacancies: {100*concentration:.1f}%")

        # Identify surface atoms (those with distance > 85% of radius)
        surface_threshold = 0.85 * self.radius
        surface_indices = []

        for i, atom in enumerate(self.atoms):
            distance = np.linalg.norm(atom.position)
            if distance >= surface_threshold:
                surface_indices.append(i)

        print(f"Identified {len(surface_indices)} surface atoms")

        # Remove fraction of surface atoms
        n_remove = int(len(surface_indices) * concentration)
        if n_remove == 0 and len(surface_indices) > 0:
            n_remove = 1

        indices_to_remove = np.random.choice(surface_indices, n_remove, replace=False)
        indices_to_remove = sorted(indices_to_remove, reverse=True)

        for idx in indices_to_remove:
            del self.atoms[idx]

        # Reassign atom IDs
        for i, atom in enumerate(self.atoms, 1):
            atom.atom_id = i

        print(f"Removed {n_remove} surface atoms")

    def create_polycrystalline(self, n_grains: int = 3):
        """
        Create polycrystalline nanoparticle with multiple grains
        Uses Voronoi tessellation with REALISTIC grain boundaries

        This version creates true grain boundaries where each grain is fully rotated.
        Some atom loss at boundaries is expected and physically realistic.
        Use LAMMPS energy minimization after generation to relax grain boundaries.

        Args:
            n_grains: Number of grains (default: 3)
        """
        if n_grains <= 1:
            return

        print(f"\nCreating REALISTIC polycrystalline structure with {n_grains} grains")

        # Check if grain count is reasonable for particle size
        avg_grain_size = self.diameter_nm / n_grains
        if avg_grain_size < 5.0:
            print(
                f"WARNING: Average grain size ({avg_grain_size:.1f} nm) is very small!"
            )
            print(f"Recommended grain sizes:")
            print(f"  < 15nm:  Use 1 grain (single crystal)")
            print(f"  15-30nm: Use 1-2 grains")
            print(f"  30-50nm: Use 2-3 grains")
            print(f"  > 50nm:  Use 3-5 grains")

        # Generate grain centers - MORE CONSERVATIVE
        grain_centers = self._generate_grain_centers_safe(n_grains)

        print(f"Grain centers positioned (safe placement)")

        # Generate random rotation for each grain
        rotations = []
        print("Generating grain orientations:")
        for i in range(n_grains):
            # Random Euler angles (ZXZ convention)
            alpha = np.random.uniform(0, 2 * np.pi)
            beta = np.random.uniform(0, np.pi)
            gamma = np.random.uniform(0, 2 * np.pi)

            # Create rotation matrix
            R = self._euler_to_rotation_matrix(alpha, beta, gamma)
            rotations.append(R)

            # Report misorientation angle
            angle = np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)) * 180 / np.pi
            print(f"  Grain {i+1}: Rotation angle = {angle:.1f}°")

        # Voronoi tessellation with PROPER boundary checking
        atoms_to_keep = []
        atoms_removed_sphere = 0
        grain_assignments = []

        print("\nApplying grain rotations with boundary checks...")

        for atom in self.atoms:
            # Find nearest grain center (Voronoi cell assignment)
            distances = [np.linalg.norm(atom.position - gc) for gc in grain_centers]
            grain_id = np.argmin(distances)

            # Apply FULL rotation relative to grain center
            relative_pos = atom.position - grain_centers[grain_id]
            rotated_pos = rotations[grain_id] @ relative_pos
            new_position = rotated_pos + grain_centers[grain_id]

            # CRITICAL FIX: Check if new position is inside the SHAPE (not just radius)
            if self._is_inside_shape(new_position):
                atom.position = new_position
                atoms_to_keep.append(atom)
                grain_assignments.append(grain_id)
            else:
                atoms_removed_sphere += 1

        # Calculate statistics
        initial_count = len(self.atoms)
        self.atoms = atoms_to_keep

        # Renumber atoms
        for i, atom in enumerate(self.atoms, 1):
            atom.atom_id = i

        # Report statistics
        atom_loss_percent = 100.0 * atoms_removed_sphere / initial_count

        print(f"\n{'='*60}")
        print("POLYCRYSTALLINE STRUCTURE STATISTICS:")
        print(f"{'='*60}")
        print(f"Number of grains:            {n_grains}")
        print(f"Average grain size:          {avg_grain_size:.1f} nm")
        print(f"Initial atoms:               {initial_count}")
        print(
            f"Atoms removed (boundaries):  {atoms_removed_sphere} ({atom_loss_percent:.1f}%)"
        )
        print(f"Final atom count:            {len(self.atoms)}")
        print(f"Retention rate:              {100-atom_loss_percent:.1f}%")

        # Count atoms per grain
        from collections import Counter

        grain_counts = Counter(grain_assignments)
        print(f"\nAtoms per grain:")
        for grain_id in range(n_grains):
            count = grain_counts.get(grain_id, 0)
            pct = 100.0 * count / len(self.atoms) if len(self.atoms) > 0 else 0
            print(f"  Grain {grain_id+1}: {count} atoms ({pct:.1f}%)")

        print(f"{'='*60}")

        # Grain imbalance check
        if len(grain_counts) > 0:
            max_grain_size = max(grain_counts.values())
            min_grain_size = min(grain_counts.values())
            imbalance_ratio = (
                max_grain_size / min_grain_size if min_grain_size > 0 else float("inf")
            )

            if imbalance_ratio > 2.5:
                print(f"\n⚠️  WARNING: Grain size imbalance ({imbalance_ratio:.1f}x)")
                print(f"   Try different random seed or reduce grain count")
            else:
                print(f"\n✓ Grain sizes well balanced (ratio: {imbalance_ratio:.1f}x)")

        if atom_loss_percent > 15.0:
            print(f"\n⚠️  High atom loss: {atom_loss_percent:.1f}%")
            print(f"   Consider fewer grains or larger particle")
        elif atom_loss_percent < 8.0:
            print(f"\n✓ Atom loss reasonable: {atom_loss_percent:.1f}%")

        print(f"\n📌 IMPORTANT: Use energy minimization in LAMMPS!")
        print(f"   minimize 1.0e-4 1.0e-6 1000 10000")

    def _generate_grain_centers_safe(self, n_grains: int) -> List[np.ndarray]:
        """
        Generate grain centers with ULTRA-CONSERVATIVE placement
        This prevents the Pac-Man effect by keeping grains well inside

        Args:
            n_grains: Number of grain centers to generate

        Returns:
            List of grain center positions
        """
        grain_centers = []

        # MUCH MORE CONSERVATIVE: Keep grains in inner 50% of radius
        # This ensures even after rotation, atoms stay within boundary
        max_radius = self.radius * 0.5  # REDUCED from 0.8 to 0.5

        # Calculate expected grain size
        volume_per_grain = (4 / 3 * np.pi * self.radius**3) / n_grains
        avg_grain_radius = (3 * volume_per_grain / (4 * np.pi)) ** (1 / 3)

        print(f"  Maximum grain center radius: {max_radius:.1f} Å (50% of particle)")
        print(f"  Expected grain size: ~{avg_grain_radius*2/10:.1f} nm")

        if avg_grain_radius * 2 > 200:  # 20nm in Angstroms
            print(f"  ⚠️  WARNING: Grain size >20nm not typical for EDM synthesis!")
            recommended = int((self.diameter_nm / 10) ** 3)
            print(f"  Recommended: Use ~{recommended} grains for 10nm grain size")

        # Minimum separation between grains (50% of grain diameter)
        min_distance = avg_grain_radius
        max_attempts = 1000

        for i in range(n_grains):
            attempts = 0
            placed = False

            while attempts < max_attempts and not placed:
                if i == 0:
                    # First grain ALWAYS at exact center
                    new_center = np.array([0.0, 0.0, 0.0])
                    grain_centers.append(new_center)
                    placed = True
                    print(f"  Grain 1: center (0, 0, 0)")
                else:
                    # Random position within inner 50%
                    theta = np.random.uniform(0, 2 * np.pi)
                    phi = np.random.uniform(0, np.pi)
                    r = np.random.uniform(0, max_radius)

                    x = r * np.sin(phi) * np.cos(theta)
                    y = r * np.sin(phi) * np.sin(theta)
                    z = r * np.cos(phi)

                    new_center = np.array([x, y, z])

                    # Check minimum distance
                    distances = [
                        np.linalg.norm(new_center - gc) for gc in grain_centers
                    ]
                    if min(distances) >= min_distance:
                        grain_centers.append(new_center)
                        placed = True
                        print(f"  Grain {i+1}: r={r:.1f}Å (attempt {attempts+1})")

                attempts += 1

            if not placed:
                # Fallback: place near center with random offset
                offset = np.random.randn(3) * avg_grain_radius * 0.3
                grain_centers.append(offset)
                print(f"  Grain {i+1}: near center (fallback)")

        return grain_centers

    def _euler_to_rotation_matrix(
        self, alpha: float, beta: float, gamma: float
    ) -> np.ndarray:
        """
        Convert Euler angles (ZXZ convention) to rotation matrix

        Args:
            alpha: First rotation angle around Z axis (radians)
            beta: Second rotation angle around X axis (radians)
            gamma: Third rotation angle around Z axis (radians)

        Returns:
            3x3 rotation matrix
        """
        # Precompute trigonometric values
        ca, sa = np.cos(alpha), np.sin(alpha)
        cb, sb = np.cos(beta), np.sin(beta)
        cg, sg = np.cos(gamma), np.sin(gamma)

        # Build rotation matrix using ZXZ convention
        R = np.array(
            [
                [ca * cg - sa * cb * sg, -ca * sg - sa * cb * cg, sa * sb],
                [sa * cg + ca * cb * sg, -sa * sg + ca * cb * cg, -ca * sb],
                [sb * sg, sb * cg, cb],
            ]
        )

        return R

    def write_xyz(self, filename: str):
        """
        Write particle to XYZ file format

        Args:
            filename: Output filename
        """
        print(f"\nWriting to {filename}...")

        with open(filename, "w") as f:
            # Header: number of atoms
            f.write(f"{len(self.atoms)}\n")

            # Comment line
            n_ni = sum(1 for atom in self.atoms if atom.element == "Ni")
            n_ti = sum(1 for atom in self.atoms if atom.element == "Ti")
            ni_pct = 100.0 * n_ni / len(self.atoms)
            ti_pct = 100.0 * n_ti / len(self.atoms)

            f.write(
                f"NiTi nanoparticle, D={self.diameter_nm}nm, "
                f"Ni{ni_pct:.1f}Ti{ti_pct:.1f}, N={len(self.atoms)}\n"
            )

            # Atom lines: element x y z
            for atom in self.atoms:
                f.write(
                    f"{atom.element} {atom.position[0]:.6f} "
                    f"{atom.position[1]:.6f} {atom.position[2]:.6f}\n"
                )

        print(f"Successfully wrote {len(self.atoms)} atoms to {filename}")

    def write_lammps_data(self, filename: str):
        """
        Write particle to LAMMPS data file format

        Args:
            filename: Output filename
        """
        print(f"\nWriting LAMMPS data file to {filename}...")

        # Assign atom types: Ni=1, Ti=2
        type_map = {"Ni": 1, "Ti": 2}

        # Calculate actual particle extent
        if len(self.atoms) > 0:
            positions = np.array([atom.position for atom in self.atoms])
            max_extent = np.max(np.abs(positions))
        else:
            max_extent = self.radius

        # Calculate box bounds with GENEROUS padding for blob shapes
        # For blob, we need extra room since it can extend to 1.25x radius
        if self.shape == "blob":
            padding = max_extent * 0.3  # 30% padding for blob
        else:
            padding = 20.0  # Standard 20 Angstrom padding

        box_size = max_extent * 2 + 2 * padding
        xlo, xhi = -box_size / 2, box_size / 2
        ylo, yhi = -box_size / 2, box_size / 2
        zlo, zhi = -box_size / 2, box_size / 2

        print(f"  Box size: {box_size:.1f} Å (particle extent: {max_extent*2:.1f} Å)")
        print(f"  Padding: {padding:.1f} Å")

        # Get composition
        n_ni = sum(1 for atom in self.atoms if atom.element == "Ni")
        n_ti = sum(1 for atom in self.atoms if atom.element == "Ti")
        ni_pct = 100.0 * n_ni / len(self.atoms) if len(self.atoms) > 0 else 0
        ti_pct = 100.0 * n_ti / len(self.atoms) if len(self.atoms) > 0 else 0

        with open(filename, "w") as f:
            # Header
            f.write(
                f"LAMMPS data file - NiTi nanoparticle D={self.diameter_nm}nm "
                f"shape={self.shape} Ni{ni_pct:.1f}Ti{ti_pct:.1f}\n\n"
            )

            # System info
            f.write(f"{len(self.atoms)} atoms\n")
            f.write("2 atom types\n\n")

            # Box bounds
            f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
            f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
            f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n\n")

            # Masses (Ni: 58.6934, Ti: 47.867 amu)
            f.write("Masses\n\n")
            f.write("1 58.6934  # Ni\n")
            f.write("2 47.8670  # Ti\n\n")

            # Atoms (atom-ID atom-type x y z)
            f.write("Atoms # atomic\n\n")
            for atom in self.atoms:
                atom_type = type_map[atom.element]
                f.write(
                    f"{atom.atom_id} {atom_type} "
                    f"{atom.position[0]:.6f} {atom.position[1]:.6f} "
                    f"{atom.position[2]:.6f}\n"
                )

        print(f"Successfully wrote LAMMPS data file with {len(self.atoms)} atoms")

    def make_amorphous(self, disorder_strength: float = 0.3):
        """
        Create amorphous structure by randomly displacing atoms
        This simulates the rapid quenching from EDM synthesis

        Args:
            disorder_strength: Amount of atomic displacement (0-1, default: 0.3)
                              0.3 means ~30% of lattice parameter displacement
        """
        print(f"\nCreating amorphous structure (rapid quench simulation)")
        print(f"Disorder strength: {disorder_strength:.2f}")

        max_displacement = self.lattice_param * disorder_strength

        for atom in self.atoms:
            # Random displacement in all directions
            displacement = np.random.randn(3) * max_displacement / 3.0
            atom.position += displacement

        print(f"Atoms randomly displaced by ~{max_displacement:.2f} Å")
        print(f"Structure is now amorphous (no long-range order)")


def main():
    """Main function with command-line interface"""
    parser = argparse.ArgumentParser(
        description="Generate spherical NiTi nanoparticles with defects",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
NiTi Nanoparticle Generator for EDM Synthesis Simulation
=========================================================

Generates realistic NiTi (Nitinol) nanoparticles for molecular dynamics simulations,
with focus on EDM (Electrical Discharge Machining) synthesis conditions.

EDM Synthesis Physics:
----------------------
- Spark plasma melts/vaporizes NiTi at 2000-3000°C
- Rapid quenching in DI water at 10^6 to 10^9 K/s
- Results in amorphous or ultra-fine nanocrystalline structures
- Typical grain sizes: 2-20 nm (not large grains!)

Available Parameters (13 total):
--------------------------------
Required:
  --diameter, -d    Nanoparticle diameter in nanometers
  --output, -o      Output filename (.xyz or .data)

Composition:
  --ni-percent      Percentage of Ni atoms (0-100, default: 50.0 for equiatomic)
  --lattice         B2 lattice parameter in Angstroms (default: 3.015)

Shape:
  --shape           Particle shape: sphere, faceted, ellipsoid, rough, or blob (default: sphere)
  --aspect-ratio    For ellipsoid shape (default: 1.0)

Structure (MUTUALLY EXCLUSIVE - use only ONE):
  --grains          Number of grains for polycrystalline (default: 1 = single crystal)
  --amorphous       Disorder strength for amorphous structure (0.0-1.0, e.g., 0.3 for EDM)

Defects:
  --vacancies       Random vacancy concentration (0.0-0.5, default: 0.0)
  --antisites       Ni-Ti swap concentration (0.0-0.5, default: 0.0)
  --surface-vac     Surface vacancy concentration (0.0-0.2, default: 0.0)

Other:
  --seed            Random seed for reproducibility
  --lammps          Also generate LAMMPS data file

Shape Options Explained (Based on SEM Analysis):
------------------------------------------------
  sphere    - Perfect mathematical sphere (baseline, academic, easy to analyze)
  ellipsoid - Elongated/stretched particle (common in flow-based synthesis)
  faceted   - Wulff-like with {100} and {110} crystal planes (slow-cooled equilibrium crystals)
  rough     - Atomic-scale surface roughness (good for interface MD studies)
  blob      - Irregular, smooth organic shape (MOST REALISTIC for EDM rapid quench!)

For EDM Synthesis: Use --shape blob to match experimental SEM images showing irregular potato-like morphology

Structure Types (MUTUALLY EXCLUSIVE - choose only ONE):
--------------------------------------------------------
  Single Crystal    - Default (no --grains or --grains 1) - Best for small particles <15nm
  Nanocrystalline   - Use --grains N where N = (diameter/10)^3 for 10nm grain size
  Amorphous         - Use --amorphous 0.3 (most realistic for EDM rapid quench)

IMPORTANT: Do NOT use --grains and --amorphous together!
  --grains     = Multiple oriented crystals (polycrystalline)
  --amorphous  = No crystal structure (glass-like)
  These are physically incompatible structure types.

Examples (All 13 Parameters Demonstrated):
------------------------------------------
# Basic single crystal sphere
python niti_np_gen.py --diameter 15 --output niti_15nm.xyz

# MOST REALISTIC EDM: Blob shape with amorphous structure
python niti_np_gen.py --diameter 60 --shape blob --amorphous 0.3 --output niti_60nm_edm_realistic.xyz --lammps

# Large EDM particle matching SEM (130nm blob, amorphous)
python niti_np_gen.py --diameter 130 --shape blob --amorphous 0.3 --vacancies 0.01 --surface-vac 0.03 --output niti_130nm_edm.xyz --lammps --seed 42

# Nanocrystalline blob (realistic internal structure)
python niti_np_gen.py --diameter 100 --shape blob --grains 125 --vacancies 0.01 --output niti_100nm_nc_blob.xyz --lammps

# Faceted single crystal (for comparison, NOT realistic for EDM)
python niti_np_gen.py --diameter 20 --shape faceted --output niti_20nm_faceted.xyz --lammps

# Ellipsoidal particle with aspect ratio and defects
python niti_np_gen.py --diameter 25 --shape ellipsoid --aspect-ratio 1.3 --vacancies 0.02 --antisites 0.01 --output niti_25nm_ellipsoid.xyz --lammps

# Ni-rich composition with all defect types
python niti_np_gen.py --diameter 50 --shape blob --grains 125 --ni-percent 52 --vacancies 0.01 --antisites 0.005 --surface-vac 0.03 --output niti_50nm_nirich.xyz --lammps

# ALL 13 PARAMETERS USED: Maximum complexity
python niti_np_gen.py --diameter 40 --shape blob --grains 64 --ni-percent 51 --lattice 3.015 --vacancies 0.01 --antisites 0.005 --surface-vac 0.03 --output niti_40nm_all_params.xyz --lammps --seed 42

# Testing different shapes (same size, same seed for comparison)
python niti_np_gen.py --diameter 30 --shape sphere --output niti_30nm_sphere.xyz --seed 100
python niti_np_gen.py --diameter 30 --shape faceted --output niti_30nm_faceted.xyz --seed 100
python niti_np_gen.py --diameter 30 --shape ellipsoid --aspect-ratio 1.5 --output niti_30nm_ellipsoid.xyz --seed 100
python niti_np_gen.py --diameter 30 --shape rough --output niti_30nm_rough.xyz --seed 100
python niti_np_gen.py --diameter 30 --shape blob --output niti_30nm_blob.xyz --seed 100

# Extreme quench simulation (highly disordered amorphous)
python niti_np_gen.py --diameter 35 --shape blob --amorphous 0.4 --vacancies 0.02 --surface-vac 0.08 --output niti_35nm_extreme_quench.xyz --lammps

Recommended Settings by Size (Based on EDM Physics):
----------------------------------------------------
  5-15nm    Single crystal blob (too small for grains)
            python niti_np_gen.py --diameter 10 --shape blob --output niti_10nm.xyz --lammps

  15-30nm   Amorphous blob (most realistic for EDM) OR single crystal
            python niti_np_gen.py --diameter 25 --shape blob --amorphous 0.3 --output niti_25nm_edm.xyz --lammps
            python niti_np_gen.py --diameter 25 --shape blob --output niti_25nm_single.xyz --lammps

  30-60nm   Amorphous blob with defects (realistic EDM conditions)
            python niti_np_gen.py --diameter 50 --shape blob --amorphous 0.3 --vacancies 0.01 --output niti_50nm_edm.xyz --lammps

  60-150nm  Nanocrystalline blob OR amorphous (matching SEM observations)
            python niti_np_gen.py --diameter 100 --shape blob --grains 125 --vacancies 0.01 --output niti_100nm_nc.xyz --lammps
            python niti_np_gen.py --diameter 130 --shape blob --amorphous 0.3 --vacancies 0.01 --output niti_130nm_amor.xyz --lammps

LAMMPS Post-Processing (CRITICAL STEP):
---------------------------------------
After generation, ALWAYS use energy minimization in LAMMPS to relax the structure:
  minimize 1.0e-4 1.0e-6 1000 10000

For polycrystalline structures, also equilibrate at finite temperature:
  fix 1 all nvt temp 300.0 300.0 0.1
  run 10000

For amorphous structures, perform annealing cycle:
  fix 1 all nvt temp 300.0 600.0 0.1
  run 5000
  fix 1 all nvt temp 600.0 300.0 0.1
  run 5000

For Shape Memory Effect (SME) and Superelasticity (SE) Studies:
---------------------------------------------------------------
1. Start with single crystal blob for baseline behavior (simplest case)
   python niti_np_gen.py --diameter 15 --shape blob --output baseline_15nm.xyz --lammps

2. Add defects progressively to study their effects on martensitic transformation
   python niti_np_gen.py --diameter 15 --shape blob --vacancies 0.01 --output with_vac_15nm.xyz --lammps
   python niti_np_gen.py --diameter 15 --shape blob --vacancies 0.01 --antisites 0.005 --output with_defects_15nm.xyz --lammps

3. Use amorphous blob to simulate REALISTIC EDM conditions (matches SEM!)
   python niti_np_gen.py --diameter 60 --shape blob --amorphous 0.3 --output edm_realistic_60nm.xyz --lammps
   python niti_np_gen.py --diameter 130 --shape blob --amorphous 0.3 --output edm_realistic_130nm.xyz --lammps

4. Compare nanocrystalline vs amorphous transformation behavior
   python niti_np_gen.py --diameter 50 --shape blob --grains 125 --output nc_50nm.xyz --lammps
   python niti_np_gen.py --diameter 50 --shape blob --amorphous 0.3 --output amor_50nm.xyz --lammps

5. Study size effects on SME/SE (use same structure, vary only size)
   python niti_np_gen.py --diameter 10 --shape blob --amorphous 0.3 --output size_10nm.xyz --lammps --seed 42
   python niti_np_gen.py --diameter 20 --shape blob --amorphous 0.3 --output size_20nm.xyz --lammps --seed 42
   python niti_np_gen.py --diameter 30 --shape blob --amorphous 0.3 --output size_30nm.xyz --lammps --seed 42
   python niti_np_gen.py --diameter 50 --shape blob --amorphous 0.3 --output size_50nm.xyz --lammps --seed 42

6. Study composition effects on transformation temperature
   python niti_np_gen.py --diameter 30 --shape blob --ni-percent 49 --amorphous 0.3 --output comp_49ni.xyz --lammps
   python niti_np_gen.py --diameter 30 --shape blob --ni-percent 50 --amorphous 0.3 --output comp_50ni.xyz --lammps
   python niti_np_gen.py --diameter 30 --shape blob --ni-percent 51 --amorphous 0.3 --output comp_51ni.xyz --lammps
   python niti_np_gen.py --diameter 30 --shape blob --ni-percent 52 --amorphous 0.3 --output comp_52ni.xyz --lammps

Note on Shape Selection:
------------------------
Based on SEM analysis of EDM-synthesized NiTi nanoparticles (60-156nm range):
- Particles show IRREGULAR, SMOOTH morphology (potato/pebble-like)
- NO crystallographic facets visible (even at 130-156nm size!)
- Surface tension from liquid droplet + rapid quench → blob shape
- Faceted shape implies slow cooling and is INCORRECT for EDM

Use --shape blob for publication-quality realistic EDM simulations!
Use --shape faceted only for academic comparison or slow-cooled reference cases.
Use --shape sphere only for idealized baseline studies.

Common Workflows:
-----------------
# Quick test (small, fast generation)
python niti_np_gen.py --diameter 10 --shape blob --output test.xyz

# Production EDM simulation (realistic, all parameters)
python niti_np_gen.py --diameter 60 --shape blob --amorphous 0.3 --ni-percent 50.5 --vacancies 0.01 --surface-vac 0.03 --output production_60nm.xyz --lammps --seed 12345

# Batch generation for size series
for size in 10 20 30 40 50; do python niti_np_gen.py --diameter $size --shape blob --amorphous 0.3 --output niti_${size}nm.xyz --lammps --seed 42; done

# Parametric study of grain count
python niti_np_gen.py --diameter 50 --shape blob --grains 27 --output grains_27.xyz --lammps
python niti_np_gen.py --diameter 50 --shape blob --grains 64 --output grains_64.xyz --lammps
python niti_np_gen.py --diameter 50 --shape blob --grains 125 --output grains_125.xyz --lammps
        """,
    )

    # Required arguments
    parser.add_argument(
        "--diameter",
        "-d",
        type=float,
        required=True,
        help="Nanoparticle diameter in nanometers",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output filename (.xyz or .data)",
    )

    # Optional arguments
    parser.add_argument(
        "--ni-percent",
        type=float,
        default=50.0,
        help="Percentage of Ni atoms (0-100, default: 50.0 for equiatomic)",
    )
    parser.add_argument(
        "--lattice",
        type=float,
        default=3.015,
        help="B2 lattice parameter in Angstroms (default: 3.015)",
    )
    parser.add_argument(
        "--vacancies",
        type=float,
        default=0.0,
        help="Vacancy concentration (0.0-0.1, default: 0.0)",
    )
    parser.add_argument(
        "--antisites",
        type=float,
        default=0.0,
        help="Antisite defect concentration (0.0-0.1, default: 0.0)",
    )
    parser.add_argument(
        "--surface-vac",
        type=float,
        default=0.0,
        help="Surface vacancy concentration (0.0-0.2, default: 0.0)",
    )
    parser.add_argument(
        "--grains",
        type=int,
        default=1,
        help="Number of grains for polycrystalline (default: 1)",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--lammps", action="store_true", help="Also generate LAMMPS data file"
    )
    parser.add_argument(
        "--shape",
        type=str,
        default="sphere",
        choices=["sphere", "ellipsoid", "faceted", "rough", "blob"],
        help="Particle shape (default: sphere)",
    )
    parser.add_argument(
        "--aspect-ratio",
        type=float,
        default=1.0,
        help="Aspect ratio for ellipsoid shape (default: 1.0)",
    )
    parser.add_argument(
        "--amorphous",
        type=float,
        default=0.0,
        help="Make amorphous structure (disorder 0.0-1.0, e.g., 0.3 for EDM)",
    )

    args = parser.parse_args()

    # Set random seed if provided
    if args.seed is not None:
        np.random.seed(args.seed)
        print(f"Random seed set to: {args.seed}")

    # Validate inputs
    if args.diameter <= 0:
        print("ERROR: Diameter must be positive")
        sys.exit(1)

    if not (0 < args.ni_percent < 100):
        print("ERROR: Ni percentage must be between 0 and 100")
        sys.exit(1)

    if not (0 <= args.vacancies < 0.5):
        print("ERROR: Vacancy concentration must be in [0, 0.5)")
        sys.exit(1)

    if not (0 <= args.antisites < 0.5):
        print("ERROR: Antisite concentration must be in [0, 0.5)")
        sys.exit(1)

    # NEW: Check for conflicting structure options
    if args.grains > 1 and args.amorphous > 0:
        print("ERROR: Cannot use --grains and --amorphous together")
        print("  --grains creates polycrystalline (multiple oriented crystals)")
        print("  --amorphous creates disordered structure (no crystal grains)")
        print("  These are mutually exclusive structure types.")
        print("\nChoose ONE:")
        print("  For nanocrystalline: --grains N (e.g., --grains 27)")
        print("  For amorphous:       --amorphous 0.3")
        sys.exit(1)

    # Create nanoparticle
    print("=" * 60)
    print("NiTi NANOPARTICLE GENERATOR")
    print("=" * 60)

    particle = NiTiNanoparticle(
        args.diameter, args.ni_percent, args.lattice, args.shape, args.aspect_ratio
    )

    # Generate base structure
    particle.generate_b2_structure()

    # Apply structure modification (ONLY ONE of these will execute)
    if args.amorphous > 0:
        # Amorphous structure (rapid quench simulation)
        particle.make_amorphous(args.amorphous)
    elif args.grains > 1:
        # Polycrystalline structure
        particle.create_polycrystalline(args.grains)
    # else: single crystal (no modification)

    # Add defects (these work with any structure type)
    if args.vacancies > 0:
        particle.add_vacancies(args.vacancies)

    if args.antisites > 0:
        particle.add_antisite_defects(args.antisites)

    if args.surface_vac > 0:
        particle.add_surface_vacancies(args.surface_vac)

    # Write output
    if args.output.endswith(".xyz"):
        particle.write_xyz(args.output)
    elif args.output.endswith(".data"):
        particle.write_lammps_data(args.output)
    else:
        # Default to XYZ
        particle.write_xyz(args.output)

    # Also write LAMMPS if requested
    if args.lammps:
        lammps_file = args.output.replace(".xyz", ".data")
        if not lammps_file.endswith(".data"):
            lammps_file += ".data"
        particle.write_lammps_data(lammps_file)

    print("\n" + "=" * 60)
    print("GENERATION COMPLETE!")
    print("=" * 60)


if __name__ == "__main__":
    main()
