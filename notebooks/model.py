import pandas as pd
from data_processing import scale_bead_chain
from initial_structures_defs import *
from openmm.app import PDBxFile, ForceField, Simulation, DCDReporter, StateDataReporter
import openmm as mm
import openmm.unit as u
from sys import stdout


class Model():
    def __init__(self, n_beads, n_connections, data_file_path = "data/ENCFF780PGS.bedpe", chrom = "chr1", sep = "\t") -> None:
        self.n_beads =  n_beads
        self.n_connections = n_connections
        self.data_file_path = data_file_path
        self.chrom = chrom
        self.sep = sep
        self.expected_loops = None


    @property
    def _get_prepared_data(self):
        df = pd.read_csv("data/ENCFF780PGS.bedpe", sep=self.sep, header= None)
        df.columns = ["chrom1", "start1", "end1", "chrom2","start2", "end2", "score"]

        df_chr1 = df[df['chrom1'] == self.chrom]
        df_chr1['middle1'] = (df_chr1['end1'] + df_chr1['start1'])/2
        df_chr1['middle2'] = (df_chr1['end2'] + df_chr1['start2'])/2

        middle_points = df_chr1[['middle1', 'middle2']].reset_index(drop=True)
        return scale_bead_chain(middle_points, self.n_connections, new_min=1, new_max=self.n_beads)
    
    def run_simulation(self):
        # 0. Generate some initial structure
        points = helisa(self.n_beads)
        write_mmcif(points,'init_struct.cif')
        generate_psf(self.n_beads,'LE_init_struct.psf')

        # 1. Define System
        pdb = PDBxFile('init_struct.cif')
        forcefield = ForceField('forcefields/classic_sm_ff.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedCutoff=1*u.nanometer)
        integrator = mm.LangevinIntegrator(310, 0.05, 100 * mm.unit.femtosecond)

        # 2. Define the forcefield
        # 2.1. Harmonic bond borce between succesive beads
        bond_force = mm.HarmonicBondForce()
        system.addForce(bond_force)
        for i in range(system.getNumParticles() - 1):
            bond_force.addBond(i, i + 1, 0.1, 3000.0)

        # Connecting selected beads
        df_scaled = self._get_prepared_data
        self.expected_loops = df_scaled

        for i in range(len(df_scaled)):
            middle1, middle2 = df_scaled.iloc[i,:]
            bond_force.addBond(middle1-1, middle2-1, 0.1, 10000)


        #2.2. Harmonic angle force between successive beads so as to make chromatin rigid
        angle_force = mm.HarmonicAngleForce()
        system.addForce(angle_force)
        for i in range(system.getNumParticles() - 2):
            angle_force.addAngle(i, i + 1, i + 2, np.pi, 0.0001)
            
        # 3. Minimize energy
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.reporters.append(StateDataReporter(stdout, 10, step=True, totalEnergy=True, potentialEnergy=True, temperature=True))
        simulation.reporters.append(DCDReporter('stochastic_LE.dcd', 10))
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(tolerance=0.001)
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open('minimized.cif', 'w')) # save minimized file

        # 4. Run md simulation
        simulation.context.setVelocitiesToTemperature(310, 0)
        simulation.step(10000)
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(pdb.topology, state.getPositions(), open('after_sim.cif', 'w')) # save minimized file

    def visualize_data(self):
        pass