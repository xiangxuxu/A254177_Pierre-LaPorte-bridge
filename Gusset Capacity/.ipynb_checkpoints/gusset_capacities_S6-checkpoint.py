# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: A254177_env
#     language: python
#     name: a254177_env
# ---

# %% [markdown]
# <img
# src="Images/COWI_logo_CMYK_orange_2000x350px-01.png" width="150"
# />
#
# ### A254177: Pierre Laporte Gusset Capacity Investigation 
# # Gusset Capacities
# **Author:**  LIXU  
# **Date:** 2024-AUG-13  
# **Checker:**    
# **Date:** 2024-XX-XX  
#
# ## Background Information
#
# In this calcbook, the gusset plate connection capacities will be calculated
# as well as the fastener resistances. These will be used in separate calcbooks
# to calculate the D/C of the gusset connections according to CSA S6-19.
#
# Gussets will be analysed for shear, tension and compression. Additonally, bolt
# resistances will be calculated.
#
#

# %% [markdown]
# ## Python Imports

# %%
import sys
import os

file_dir = os.path.abspath('')
parent_dir = os.path.dirname(file_dir)
sys.path.append(parent_dir)

from gusset_funcs_S6 import (in_to_m, in_to_mm,
                          E, G_s,
                          f_u_bolt)


import englib.csa.s6_19 as s6
import pandas as pd
import numpy as np



# %% [markdown]
# ## Section 14 Factors
#
# First, toggle to include Section 14 resistance adjustment factors or not. These
# should not be used to determine the acceptable section loss. For more
# information, refer to CSA S6-19, Cl. 14.14.2.
#
# <img
# src="Images/table_14_10.png" width="600"
# />
#
# Set `include_U` to `True` is you want to include the factors and to `False`, if
# not.

# %%
include_U = False

if include_U==True:
    U_net_ten = 1.18
    U_gross = 1.01
    U_bolt = 1.20
    output_suffix = "_with_resistance_adjustment.csv"
else:
    U_net_ten = 1.00
    U_gross = 1.00
    U_bolt = 1.00
    output_suffix = ".csv"


# %% [markdown]
# ## Define Gussets

# %%
# Get the current script's directory
current_dir = os.getcwd()

# Construct the full path to the CSV file
csv_file_path = os.path.join(current_dir, 'gusset_input.csv')

# Read the CSV file into a DataFrame
# Use the `index_col` parameter to specify the column(s) to set as the index
gus_df = pd.read_csv(csv_file_path,  index_col=0)

gus_df

# %%


# 'gus_df' has an index 
index_names = gus_df.index.astype(str)

# Determine n_gussets =1 for gussets connect with top bracings, and n_gussets=2 for other connections
n_gussets = np.where(index_names.str.contains('TB'), 1, 2)

# Determine Fu=450 MPa for material grade is G40.8 and 480 MPa for material grade is A441 
f_u_exist = np.where(index_names.str.contains('G40.8'), 450, 480) # in MPa 

# Determine Fy=280 MPa for material grade is G40.8 and 350 MPa for material grade is A441 
f_y_exist = np.where(index_names.str.contains('G40.8'), 280, 350) # in MPa 


# %% [markdown]
# ## Fastener Resistances
#
# ### Bolt Bearing
# Note, governed by the thinnest member flange.

# %%
# bolt bearing resistance of diagonal member 
B_r_bolt_diag = s6.equation_10_18_2_3_3_a(φ_br=U_bolt*0.8,
                                     n=np.minimum(gus_df['diag_1_bolts'],gus_df['diag_2_bolts'])*n_gussets,
                                     d=gus_df['bolt_size'],
                                     t=np.minimum(gus_df['t_diag'],gus_df['thickness']),
                                     F_u=f_u_exist)['B_r']


gussets_capacity = pd.DataFrame(index=gus_df.index)
gussets_capacity['B_r_bolt_diag']=B_r_bolt_diag

# bolt bearing resistance of vertical member 
B_r_bolt_vert = s6.equation_10_18_2_3_3_a(φ_br=U_bolt*0.8,
                                     n=gus_df['vert_bolts']*n_gussets,
                                     d=gus_df['bolt_size'],
                                     t=np.minimum(gus_df['t_vert'],gus_df['thickness']),
                                     F_u=f_u_exist)['B_r']
gussets_capacity['B_r_bolt_vert']=B_r_bolt_vert

gussets_capacity


# %% [markdown]
# ### Bolt Shear
# Conservatively assume that threads are intercepted by the shear plane. Note
# that both top and btm gussets are 5/8" thick, will only calculate once.
#

# %%
# bolt shear resistance of diagonal member 
V_r_bolt_diag = s6.equation_10_18_2_3_3_b(φ_b=U_bolt*0.8,
                                     n=np.minimum(gus_df['diag_1_bolts'],gus_df['diag_2_bolts'])*n_gussets,
                                     m=1,
                                     A_b=np.pi*gus_df['bolt_size']**2/4,
                                     F_u=f_u_bolt,
                                     longsplice=False,
                                     intercepted=True)['V_r']

gussets_capacity['V_r_bolt_diag'] = V_r_bolt_diag

# bolt shear resistance of vertical member 
V_r_bolt_vert = s6.equation_10_18_2_3_3_b(φ_b=U_bolt*0.8,
                                     n=gus_df['vert_bolts']*n_gussets,
                                     m=1,
                                     A_b=np.pi*gus_df['bolt_size']**2/4,
                                     F_u=f_u_bolt,
                                     longsplice=False,
                                     intercepted=True)['V_r']

gussets_capacity['V_r_bolt_vert'] = V_r_bolt_vert



gussets_capacity



# %% [markdown]
# ## Gusset Shear
#
# Check full depth of top connections using the resistance factors presented in S6-25

# %%

# CSA S6-25 resistance factors
Φ_gv = 0.8

# %% [markdown]
#
# ### Shear Yielding - Full Plane

# %%
# Calculate horizontal gross shear area
A_vg_h = gus_df['width']*gus_df['thickness']

# Calculate vertical gross shear area
A_vg_v = gus_df['height']*gus_df['thickness']
        
# check horizontal section resistance
gussets_capacity['V_ry_h'] = Φ_gv*0.5*f_y_exist*A_vg_h*n_gussets

# check vertical section resistance
gussets_capacity['V_ry_v'] = Φ_gv*0.5*f_y_exist*A_vg_v*n_gussets



# %% [markdown]
# ### Shear Rupture

# %%
# Calculate A_vn_h and A_vn_v for each instance in type1_instances
A_vn_h = A_vg_h - gus_df['net_sect_bolts_h']*gus_df['bolt_hole_size']*gus_df['thickness']
A_vn_v = A_vg_v - gus_df['net_sect_bolts_v']*gus_df['bolt_hole_size']*gus_df['thickness']

# check horizontal section resistance

gussets_capacity['V_ru_h'] = Φ_gv*0.5*f_u_exist* A_vn_h*n_gussets

# check vertical section resistance

gussets_capacity['V_ru_v'] = Φ_gv*0.5*f_u_exist*A_vn_v*n_gussets
gussets_capacity


# %% [markdown]
# ## Gusset Compression
#
# Calculate the compression capacity of the gusset based on the Whitmore area.
#
#
#
# <img
# src="Images/whitmore.png" width="600"
# />

# %%
φ_gc = 0.75  # resistance factor for gusset compression from CSA S6-25

# %%
#Compression capacity of the gusset plate for the diagonal member
gus_df['diag_A_whit'] = gus_df['diag_w_whit']*gus_df['thickness']

gussets_capacity['C_r_diag'] = n_gussets*s6.equation_10_9_3_1(E_s=E,
                                    r=gus_df['thickness']/(12**0.5),
                                    A=gus_df['diag_A_whit'],
                                    K=0.5, # based on MBE and CSA S6-25
                                    L=gus_df['diag_L_mid'],
                                    n=1.34,
                                    φ_s=φ_gc,
                                    F_y=f_y_exist)['C_r']

#Compression capacity of the gusset plate for the vertical member
gus_df['vert_A_whit'] = gus_df['vert_w_whit']*gus_df['thickness']

gussets_capacity['C_r_vert'] = n_gussets*s6.equation_10_9_3_1(E_s=E,
                                    r=gus_df['thickness']/(12**0.5),
                                    A=gus_df['vert_A_whit'],
                                    K=0.5, # based on MBE and CSA S6-25
                                    L=gus_df['vert_L_mid'],
                                    n=1.34,
                                    φ_s=φ_gc,
                                    F_y=f_y_exist)['C_r']

# %% [markdown]
# ## Gusset Tension
#
# Going to check the tension capacity of the gusset across the whitmore section,
# using the equations from CSA S6-19
#
# Note, using CSA equations here and therefore using the corresponding CSA
# resistance factors

# %%
# tensile yielding capacity of gusset plate for diagonal member CSA S6-19.10.8.2 (a)

T_ry_diag = n_gussets*s6.equation_10_8_2_a(φ_s=U_gross*0.95,
                                      A_g=gus_df['diag_A_whit'],
                                      F_y=f_y_exist)['T_r']

# tensile fracutre capacity of gusset plate for diagonal member CSA S6-19.10.8.2 (b)
diag_A_n = (gus_df['diag_w_whit']-gus_df['diag_net_sect_bolts_whit']*gus_df['bolt_hole_size'])*gus_df['thickness']
T_ru_diag = n_gussets*s6.equation_10_8_2_c(φ_u=U_net_ten*0.8,
                                     A_ne=diag_A_n,
                                     F_u=f_u_exist)['T_r']


# block shear capacity of gusset plate for diagonal member CSA S6-19.10.8.1.3.2.5                               
diag_A_n_bs = (gus_df['diag_width_bs'] - gus_df['diag_net_sect_bolts_bs'] * gus_df['bolt_hole_size'])*gus_df['thickness']
diag_A_gv = 2*gus_df['diag_length_bs']*gus_df['thickness']

T_r_bs_diag = n_gussets*s6.equation_10_8_1_3_2_5(φ_u=0.8,
                                           U_t=1.0,
                                           A_n=diag_A_n_bs,
                                           A_gv=diag_A_gv,
                                           F_u=f_u_exist,
                                           F_y=f_y_exist)['T_r']


# find governing tension capacity for diagonal member 
gussets_capacity['T_r_diag'] = np.minimum(np.minimum(T_ry_diag, T_ru_diag), T_r_bs_diag)

# tensile yielding capacity of gusset plate for vertical member CSA S6-19.10.8.2 (a)

T_ry_vert = n_gussets*s6.equation_10_8_2_a(φ_s=U_gross*0.95,
                                      A_g=gus_df['vert_A_whit'],
                                      F_y=f_y_exist)['T_r']

# tensile fracutre capacity of gusset plate for vertical member CSA S6-19.10.8.2 (b)
vert_A_n = (gus_df['vert_w_whit']-gus_df['vert_net_sect_bolts_whit']*gus_df['bolt_hole_size'])*gus_df['thickness']
T_ru_vert = n_gussets*s6.equation_10_8_2_c(φ_u=U_net_ten*0.8,
                                     A_ne=vert_A_n,
                                     F_u=f_u_exist)['T_r']


# block shear capacity of gusset plate for vertical member CSA S6-19.10.8.1.3.2.5                               
vert_A_n_bs = (gus_df['vert_width_bs'] - gus_df['vert_net_sect_bolts_bs'] * gus_df['bolt_hole_size'])*gus_df['thickness']
vert_A_gv = 2*gus_df['vert_length_bs']*gus_df['thickness']

T_r_bs_vert = n_gussets*s6.equation_10_8_1_3_2_5(φ_u=0.8,
                                           U_t=1.0,
                                           A_n=vert_A_n_bs,
                                           A_gv=vert_A_gv,
                                           F_u=f_u_exist,
                                           F_y=f_y_exist)['T_r']


# find governing tension capacity for vertical member 
gussets_capacity['T_r_vert'] = np.minimum(np.minimum(T_ry_vert, T_ru_vert), T_r_bs_vert)
gussets_capacity =gussets_capacity/1000


#Round the values to 1 decimal place
gussets_capacity = gussets_capacity.round(1)

# reorder columns
desired_order = ['B_r_bolt_diag', 'V_r_bolt_diag','C_r_diag','T_r_diag', 'B_r_bolt_vert','V_r_bolt_vert','C_r_vert','T_r_vert','V_ry_h', 'V_ry_v', 'V_ru_h','V_ru_v']
gussets_capacity = gussets_capacity[desired_order]

gussets_capacity


# %%
# !jupyter nbconvert --to html_embed gusset_capacities_copy.ipynb

# %%

# %%

# %%

# %%

# %%
