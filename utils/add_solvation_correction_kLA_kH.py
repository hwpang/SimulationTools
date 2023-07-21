from rmgpy import settings
from rmgpy.yml import write_yml
from rmgpy.chemkin import save_chemkin, load_chemkin_file
from rmgpy.rmg.main import RMG
from rmgpy.data.thermo import find_cp0_and_cpinf
from rmgpy.thermo.thermoengine import process_thermo_data
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.data.vaporLiquidMassTransfer import (
    vapor_liquid_mass_transfer,
    liquidVolumetricMassTransferCoefficientPowerLaw,
)

species_dictionary_path = "species_dictionary.txt"
chemkin_path = "chem_annotated.inp"
save_path = "chem_liquid_kLA_kH.rms"

liqspcs, liqrxns = load_chemkin_file(
    chemkin_path,
    dictionary_path=species_dictionary_path,
    check_duplicates=False,
)


rmg = RMG()

rmg.database_directory = settings["database.directory"]
rmg.thermo_libraries = ["primaryThermoLibrary"]
rmg.kinetics_families = "default"
rmg.kinetics_depositories = ["training"]
rmg.kinetics_estimator = "rate rules"
rmg.solvent = "undecane"
rmg.reaction_libraries = []

rmg.load_database()

solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
diffusion_limiter.enable(solvent_data, rmg.database.solvation)

liquid_volumetric_mass_transfer_coefficient_power_law = (
    liquidVolumetricMassTransferCoefficientPowerLaw(
        prefactor=14.5,
        diffusion_coefficient_power=1 / 2,
        solvent_viscosity_power=-1 / 6,
        solvent_density_power=-1 / 6,
    )
)

vapor_liquid_mass_transfer.enable(
    solvent_data,
    rmg.database.solvation,
    liquid_volumetric_mass_transfer_coefficient_power_law,
)

solvliqspcs = liqspcs
for spc in solvliqspcs:
    find_cp0_and_cpinf(spc, spc.thermo)
    spc.thermo = process_thermo_data(spc, spc.thermo, solvent_name=rmg.solvent)
    spc.get_liquid_volumetric_mass_transfer_coefficient_data()
    spc.get_henry_law_constant_data()

write_yml(
    solvliqspcs,
    liqrxns,
    solvent=rmg.solvent,
    solvent_data=solvent_data,
    path=save_path,
)
