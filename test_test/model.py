import pandas as pd
from data_processing import scale_bead_chain
from initial_structures_defs import *
from openmm.app import PDBxFile, ForceField, Simulation, DCDReporter, StateDataReporter
import openmm as mm
import openmm.unit as u
from sys import stdout


class Model():
    def __init__(self, n_beads, n_connections, data_file_path = "data/ENCFF780PGS.bedpe", chrom = "chr1", sep = "\t", struct = hilbert_curve3d) -> None:
        self.n_beads =  n_beads
        self.n_connections = n_connections
        self.data_file_path = data_file_path
        self.chrom = chrom
        self.sep = sep
        self.expected_loops = None
        self.struct = struct

        self.forcefield = None
        self.system = None
        self.integrator = None
        self.pdb = None

        self.ENV_POWER = 6
        self.ENV_SIGMA = 0.1
        self.ENV_R_SMALL = 0.02
        self.ENV_EPSILON = 100
        


    @property
    def _get_prepared_data(self):
        df = pd.read_csv("data/ENCFF780PGS.bedpe", sep=self.sep, header= None)
        df.columns = ["chrom1", "start1", "end1", "chrom2","start2", "end2", "score"]

        df_chr1 = df[df['chrom1'] == self.chrom]
        df_chr1['middle1'] = (df_chr1['end1'] + df_chr1['start1'])/2
        df_chr1['middle2'] = (df_chr1['end2'] + df_chr1['start2'])/2

        middle_points = df_chr1[['middle1', 'middle2']].reset_index(drop=True)
        return scale_bead_chain(middle_points, self.n_connections, new_min=1, new_max=self.n_beads)
    
    
    def _define_init_struct(self, struct = hilbert_curve3d):
        points = struct(self.n_beads)
        write_mmcif(points,f'initial_structures_tests/{struct.__name__}/init_struct.cif')
        generate_psf(self.n_beads,f'initial_structures_tests/{struct.__name__}/LE_init_struct.psf')

    def _define_system(self):
        self.pdb = PDBxFile(f'initial_structures_tests/{self.struct.__name__}/init_struct.cif')
        self.forcefield = ForceField('forcefields/classic_sm_ff.xml')
        self.system = self.forcefield.createSystem(self.pdb.topology, nonbondedCutoff=1*u.nanometer)
        self.integrator = mm.LangevinIntegrator(310, 0.05, 100 * mm.unit.femtosecond)

    def _add_bond_force(self):
        bond_force = mm.HarmonicBondForce()
        self.system.addForce(bond_force)
        for i in range(self.n_beads - 1): # make a line 
            bond_force.addBond(i, i + 1, 0.1, 3000.0)

        # connecting selected beads
        df_scaled = self._get_prepared_data
        self.expected_loops = df_scaled

        for i in range(len(df_scaled)):
            middle1, middle2 = df_scaled.iloc[i,:]
            bond_force.addBond(middle1-1, middle2-1, 0.1, 10000)

    def _add_lennard_jonnes_force(self):
        env_force = mm.CustomNonbondedForce(f'epsilon*(sigma/(r+r_small))^{self.ENV_POWER}')
        env_force.setForceGroup(1)
        env_force.addGlobalParameter('epsilon', defaultValue=self.ENV_EPSILON)
        env_force.addGlobalParameter('r_small', defaultValue=self.ENV_R_SMALL)
        env_force.addGlobalParameter('sigma', defaultValue=self.ENV_SIGMA)

        for i in range(self.n_beads):
            env_force.addParticle()
        self.system.addForce(env_force)
    
    def _add_angle_forces(self):
        angle_force = mm.HarmonicAngleForce()
        self.system.addForce(angle_force)
        for i in range(self.n_beads - 2):
            angle_force.addAngle(i, i + 1, i + 2, np.pi, 0.0001)

    def run_simulation(self):
        self._define_init_struct(struct= self.struct)

        self._define_system()

        # Define the forcefield
        self._add_bond_force()

        self._add_lennard_jonnes_force()

        self._add_angle_forces()
            
        # 3. Minimize energy
        simulation = Simulation(self.pdb.topology, self.system, self.integrator)
        simulation.reporters.append(StateDataReporter(stdout, 10, step=True, totalEnergy=True, potentialEnergy=True, temperature=True))
        simulation.reporters.append(DCDReporter(f'initial_structures_tests/{self.struct.__name__}/stochastic_LE.dcd', 10))
        simulation.context.setPositions(self.pdb.positions)
        simulation.minimizeEnergy(tolerance=0.001)
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(self.pdb.topology, state.getPositions(), open(f'initial_structures_tests/{self.struct.__name__}/minimized.cif', 'w')) # save minimized file

        # 4. Run md simulation
        simulation.context.setVelocitiesToTemperature(310, 0)
        simulation.step(10000)
        state = simulation.context.getState(getPositions=True)
        PDBxFile.writeFile(self.pdb.topology, state.getPositions(),  open(f'initial_structures_tests/{self.struct.__name__}/after_sim.cif', 'w')) # save minimized file
        
    def visualize_data(self):
        pass