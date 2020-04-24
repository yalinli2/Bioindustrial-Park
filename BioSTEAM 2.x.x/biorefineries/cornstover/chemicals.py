# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:04 2019

@author: yoelr
"""
import thermosteam as tmo
from biorefineries.lipidcane.chemicals import lipidcane_chemicals
from thermosteam import functional as fn
import pandas as pd

__all__ = ('cornstover_chemicals', 'get_grouped_chemicals', 'chemical_groups')

# %%  Constants

# Common structural carbohydrates properties

# Assume heat capacity of lignin, cellulose, and hemicellulose
# and all components at 350 K are about the same [2,2].
Cp_cellulosic = 1.364

# Assume density is similar for most solids
rho = 1540 # kg/m3

cal2joule = 4.184

# %% Initialize Chemicals object and define functions

chems = cornstover_chemicals = tmo.Chemicals([])

# Chemicals to define
chemical_IDs = [
        'Water', 'Ethanol', 'Glucose', 'Galactose',
        'Mannose', 'Xylose', 'Arabinose', 'Cellobiose',
        'Sucrose', 'GlucoseOligomer', 'GalactoseOligomer',
        'MannoseOligomer', 'XyloseOligomer', 'ArabinoseOligomer',
        'Extract','SolubleLignin','HMF', 'Furfural', 'AceticAcid',
        'LacticAcid', 'Xylitol', 'Glycerol', 'SuccinicAcid',
        'NH3', 'H2SO4', 'NH4SO4', 'AmmoniumAcetate', 'DAP',
        'HNO3', 'NaNO3', 'NaOH', 'CellulaseNutrients',
        'Denaturant', 'Oil', 'Cellulose', 'Galactan', 'Mannan', 'Glucan',
        'Xylan', 'Arabinan', 'Lignin', 'Acetate', 'Protein',
        'Ash', 'Enzyme', 'DenaturedEnzyme', 'Z_mobilis', 'T_reesei',
        'Biomass', 'Tar', 'CaO', 'CaSO4', 'Graphite', 'N2', 'O2', 'CO2',
        'CH4', 'H2S', 'SO2', 'NO', 'CO', 'AmmoniumSulfate', 'NO2', 'CSL',
        'WWTsludge', 'Cellulase'
]

def append_single_phase_chemical(ID, search_ID=None, **data):
    chemical = tmo.Chemical(ID, search_ID=search_ID, **data)
    try: chemical.at_state(phase=chemical.phase_ref)
    except: pass
    chemical.default()    
    chems.append(chemical)

def extend_single_phase_chemicals(IDs):
    for ID in IDs: append_single_phase_chemical(ID)

def append_new_single_phase_chemical(ID, source=None, **data):
    chemical = tmo.Chemical.blank(ID, **data)
    if source: 
        default_phase_ref = source.phase_ref
        chemical.copy_models_from(source)
    else:
        default_phase_ref = 'l'
    if not chemical.phase_ref:
        chemical.phase_ref = default_phase_ref
    chemical.at_state(chemical.phase_ref)
    chemical.default()
    chems.append(chemical)

def append_chemical_copy(ID, chemical):
    new_chemical = chemical.copy(ID)
    chems.append(new_chemical)

def set_Cp(single_phase_chemical, Cp):
    chem = single_phase_chemical
    chem.Cn.add_model(Cp * chem.MW, top_priority=True)

def set_rho(single_phase_chemical, rho):
    V = fn.rho_to_V(rho, single_phase_chemical.MW)
    single_phase_chemical.V.add_model(V, top_priority=True)

# %% Define species

# TODO: Add heats of combustion

# As is in data bank
chems.extend(
    tmo.Chemicals(['Water', 'Ethanol', 'AceticAcid',
                   'Furfural', 'Glycerol', 'H2SO4',
                   'LacticAcid', 'SuccinicAcid', 'P4O10'])
)
chems.P4O10.at_state('s')
chems.H2SO4.at_state('l')
append_single_phase_chemical('HNO3', 'NitricAcid')
append_single_phase_chemical('Denaturant', 'Octane')
append_single_phase_chemical('DAP', 'Diammonium Phosphate')
append_single_phase_chemical('AmmoniumAcetate')
append_single_phase_chemical('AmmoniumSulfate')
append_single_phase_chemical('NaNO3', 'SodiumNitrate')
append_single_phase_chemical('Oil', 'Oleic acid')
append_single_phase_chemical('HMF')

# Will remain in the vapor phase
extend_single_phase_chemicals(['N2', 'NH3', 'O2', 'CH4', 'H2S', 'SO2'])
append_single_phase_chemical('CO2')

# Analagous vapors
append_new_single_phase_chemical('NO2', chems.N2,
                                 formula='NO2',
                                 Hf=7925*cal2joule)
append_new_single_phase_chemical('NO', chems.N2,
                                 formula='NO',
                                 Hf=82.05)
append_single_phase_chemical('CO', 'Carbon monoxide', Hf=-110.522)

# Will remain as  solid
extend_single_phase_chemicals(['Glucose', 'Xylose', 'Sucrose'])
append_single_phase_chemical('CaSO4')

subgroup = chems['Glucose', 'Xylose', 'Sucrose', 'CaSO4', 'AmmoniumSulfate']
for chemical in subgroup: set_Cp(chemical, Cp_cellulosic)

# Analagous sugars
append_chemical_copy('Mannose', chems.Glucose)
append_chemical_copy('Galactose', chems.Glucose)
append_chemical_copy('Arabinose', chems.Xylose)

# Other analogues
append_chemical_copy('CellulaseNutrients', chems.Glucose)
append_chemical_copy('Extract', chems.Glucose)
append_chemical_copy('Acetate', chems.AceticAcid)
append_chemical_copy('Tar', chems.Xylose)
chems.Acetate.Hf = -103373

# Chemicals taken from previous study
chems.append(lipidcane_chemicals.CaO)
chems.append(lipidcane_chemicals.Ash)
chems.append(lipidcane_chemicals.NaOH)
append_new_single_phase_chemical('Lignin',
                                 formula='C8H8O3',
                                 Hf=-108248*cal2joule)
set_rho(chems.Lignin, 1540)
set_Cp(chems.Lignin, Cp_cellulosic)
append_chemical_copy('SolubleLignin', chems.Lignin)

# Create structural carbohydrates
append_chemical_copy('GlucoseOligomer', chems.Glucose)
set_Cp(chems.GlucoseOligomer, Cp_cellulosic)
chems.GlucoseOligomer.formula = "C6H10O5"
chems.GlucoseOligomer.Hf = -233200*cal2joule

append_chemical_copy('GalactoseOligomer', chems.GlucoseOligomer)
append_chemical_copy('MannoseOligomer', chems.GlucoseOligomer)
append_chemical_copy('XyloseOligomer', chems.Xylose)
set_Cp(chems.XyloseOligomer, Cp_cellulosic)
chems.XyloseOligomer.formula = "C5H8O4"
chems.XyloseOligomer.Hf = -182100*cal2joule

append_chemical_copy('ArabinoseOligomer', chems.XyloseOligomer)

# Other
append_new_single_phase_chemical('Z_mobilis', formula="CH1.8O0.5N0.2",
                                 Hf=-31169.39*cal2joule)
append_new_single_phase_chemical('T_reesei', formula="CH1.645O0.445N0.205S0.005",
                                 Hf=-23200.01*cal2joule)
append_new_single_phase_chemical('Biomass', formula="CH1.64O0.39N0.23S0.0035",
                                 Hf=-23200.01*cal2joule)
append_new_single_phase_chemical('Cellulose', formula="C6H10O5", # Glucose monomer minus water
                                 Hf=-233200.06*cal2joule)
append_new_single_phase_chemical('Protein', formula="CH1.57O0.31N0.29S0.007",
                                 Hf=-17618*cal2joule)
append_new_single_phase_chemical('Enzyme', formula="CH1.59O0.42N0.24S0.01",
                                 Hf=-17618*cal2joule)
append_new_single_phase_chemical('Glucan', formula='C6H10O5',
                                 Hf=-233200*cal2joule)
append_new_single_phase_chemical('Xylan', formula="C5H8O4",
                                 Hf=-182100*cal2joule)
append_new_single_phase_chemical('Xylitol', formula="C5H12O5",
                                 Hf=-243145*cal2joule)
append_new_single_phase_chemical('Cellobiose', formula="C12H22O11",
                                 Hf=-480900*cal2joule)
append_new_single_phase_chemical('CSL', MW=1,
                                 Hf=(chems.Protein.Hf/4
                                     + chems.Water.Hf/2
                                     + chems.LacticAcid.Hf/4))
append_chemical_copy('DenaturedEnzyme', chems.Enzyme)
append_chemical_copy('Arabinan', chems.Xylan)
append_chemical_copy('Mannan',   chems.Glucan)
append_chemical_copy('Galactan', chems.Glucan)

# %% TODO: Maybe remove this

# For waste water
append_chemical_copy('WWTsludge', chems.Biomass)
append_chemical_copy('Cellulase', chems.Enzyme)

# %% Grouped chemicals

chemical_groups = dict(
    OtherSugars = ('Arabinose',
                   'Mannose',
                   'Galactose',
                   'Cellobiose',
                   'Sucrose'),
    SugarOligomers = ('GlucoseOligomer',
                      'XyloseOligomer',
                      'GalactoseOligomer',
                      'ArabinoseOligomer',
                      'MannoseOligomer'),
    OrganicSolubleSolids = ('AmmoniumAcetate',
                            'SolubleLignin',
                            'Extract', 
                            'LacticAcid', 
                            'Cellulase'),
    InorganicSolubleSolids = ('AmmoniumSulfate',
                              'DAP',
                              'NaOH',
                              'HNO3',
                              'NaNO3'),
    Furfurals = ('Furfural',
                 'HMF'),
    OtherOrganics = ('Glycerol',
                     'Denaturant',
                     'Oil',
                     'SuccinicAcid',
                     'Xylitol'),
    COxSOxNOxH2S = ('NO',
                    'NO2',
                    'SO2',
                    'CO',
                    'H2S'),
    Protein = ('Protein',
               'Enzyme',
               'DenaturedEnzyme'),
    CellMass = ('WWTsludge',
                'Z_mobilis',
                'T_reesei'),
    OtherInsolubleSolids = ('Tar',
                            'Ash',
                            'Lime'),
    OtherStructuralCarbohydrates = ('Arabinan', 
                                    'Mannan', 
                                    'Galactan')
)

tmo.settings.set_thermo(chems)
thermo = tmo.settings.get_thermo()
chems.set_synonym('CaO', 'Lime')
chems.set_synonym('Water', 'H2O')
chems.set_synonym('H2SO4', 'SulfuricAcid')
chems.set_synonym('NH3', 'Ammonia')
chems.set_synonym('AmmoniumSulfate', 'NH4SO4')
chems.set_synonym('Denaturant', 'Octane')
chems.set_synonym('CO2', 'CarbonDioxide')

def get_grouped_chemicals(stream, units='kmol/hr'):
    new_stream = tmo.Stream(thermo=thermo)
    new_stream.set_flow(stream.mol, units, stream.chemicals.IDs)
    data = {group: new_stream.get_flow(units, IDs).sum() for group, IDs in chemical_groups.items()}
    return pd.Series(data)

# Fix sugar properties


# cellulosic = ('Glucan', 'Xylan', 'GlucoseOligomer','XyloseOligomer',
#               'Acetate', 'SolubleLignin')

