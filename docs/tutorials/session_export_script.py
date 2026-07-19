# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     kernelspec:
#       display_name: Python 3
#       name: python3
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
# ---

# %% [markdown]
# # Session JSON and export script
#
# Both GUIs share `McfSession`.  Saving state captures paths, stretches, colors,
# and compose settings (not pixel data).  `export_script()` writes a standalone
# recreation script, including reproject/align steps when they were used.
# See [session_model](../guide/session_model.md).

# %%
import multicolorfits as mcf

# %% [markdown]
# ## Save / load

# %%
# s = mcf.McfSession()
# s.load_files(['a.fits', 'b.fits'], colors=['#f00', '#0ff'], labels=['A', 'B'])
# s.compose.combine_mode = 'lab'
# s.compose.show_combo_swatch = True
# s.save_state('my_session.json')
#
# s2 = mcf.McfSession()
# warnings = s2.load_state('my_session.json')
# assert warnings == []
#
# # Launch GUI from state:
# # mcf.gui(state='my_session.json')
# # mcf.gui_qt(state='my_session.json')
# # CLI: mcf-web --state my_session.json

# %% [markdown]
# ## Growing panel count
#
# Sessions saved with more than four panels restore full tab counts (up to 16).
# In scripts:
#
# ```python
# s.add_panel(color='#88FF88')
# s.set_n_panels(6)
# s.remove_panel(5)
# ```

# %% [markdown]
# ## Export script
#
# ```python
# text = s.export_script()
# open('recreate_plot.py', 'w').write(text)
# ```
#
# The emitted script uses public helpers (`load_fits`, stretches, compose flags,
# optional `reproject_image` replay) so coauthors can regenerate the figure without
# opening the GUI.

# %%
print('McfSession.save_state / load_state / export_script')
