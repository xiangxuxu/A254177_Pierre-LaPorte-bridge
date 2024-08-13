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
# ### A254177: Pierre Laporte Gusset Capacity Investigation 
# # Gusset Funcs
# **Author:**  LIXU  
# **Date:** 2024-AUG-13  
# ## Background Information
#
# In this script user functions and common materials are defined for the capacity
# calcualtions of the stiffening truss gusset connections.

# %%
import sys
import os

current_dir = os.getcwd()
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)


import pandas as pd
import numpy as np

# %%
# Material Props ---------------------------------------------------------------

# Existing Steel properties 


E = 200000  # in MPa
G_s = 77000  # in MPa

# Existing bolt properties
f_u_bolt = 827  # in MPa for A325 Bolts


# %%
# USER FUNCTIONS ---------------------------------------------------------------

def in_to_m(x):
    """Converts distance from inches to meters

    Args:
        x (float): distance in inches

    Returns:
        y (float): distance in meters
    """

    y = x*0.0254

    return y


# %%
# USER FUNCTIONS ---------------------------------------------------------------
def in_to_mm(x):
    """Converts distance from inches to millimeters

    Args:
        x (float): distance in inches

    Returns:
        y (float): distance in millimeters
    """

    y = x*25.4

    return y

# %%

# %%

# %%

# %%



# %%

# %%

# %%
