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
# # Embed the multicolorfits web GUI in a notebook
#
# This tutorial shows how to run the full browser GUI **inside** Jupyter,
# Google Colab, VS Code notebooks, or Binder — useful for workshops and demos.
#
# **Requirements:** `pip install multicolorfits[web]` (FastAPI + uvicorn).
#
# See also: `mcf.gui()` for a blocking desktop-style launch that opens a
# separate browser tab.

# %%
import multicolorfits as mcf

# %% [markdown]
# ## Quick embed (auto-detect environment)
#
# `gui_embed()` starts the server in a background thread and shows an inline
# iframe (Jupyter / VS Code), Colab's port proxy, or prints the URL when run
# as a plain script.

# %%
# Uncomment when FITS paths are available on this machine:
# s = mcf.McfSession()
# s.load_files(
#     ['halpha.fits', 'oiii.fits'],
#     colors=['#C11B17', '#4CC417'],
#     labels=['Hα', '[OIII]'],
# )
# mcf.gui_embed(session=s, height=950)

# %% [markdown]
# ## Google Colab
#
# Upload FITS files to the VM (or mount Drive), then force Colab display mode:
#
# ```python
# from google.colab import files
# uploaded = files.upload()  # pick local FITS files
# paths = list(uploaded.keys())
# mcf.gui_embed(files=paths, display='colab', height=900, print_url=True)
# ```
#
# **Note:** Colab paths live only in the notebook VM — preload with `files=`
# rather than paths from your laptop.

# %% [markdown]
# ## Local Jupyter / VS Code notebooks
#
# ```python
# mcf.gui_embed(session=s, display='iframe', height=900)
# ```
#
# The iframe loads `http://127.0.0.1:<port>` on the same machine as the kernel.

# %% [markdown]
# ## URL only / scripting
#
# ```python
# url = mcf.gui_embed(display='none', return_url=True)
# # or
# srv = mcf.start_gui_server(session=s)
# print(srv.url)
# ```
#
# Use `display='none'` when you want a handle to `GuiServer` without rendering,
# or `return_url=True` for automation and tests.

# %% [markdown]
# ## Binder / JupyterHub
#
# When `JUPYTERHUB_SERVICE_PREFIX` is set, `gui_embed()` automatically uses the
# hub's `/proxy/{port}/` URL in the iframe.  Override with `proxy_url=` if your
# deployment uses a custom layout.

# %% [markdown]
# ## Caveats
#
# | Topic | Guidance |
# |-------|----------|
# | **Blocking `mcf.gui()`** | Opens a tab and blocks the cell — use `gui_embed()` instead. |
# | **Port clashes** | The launcher picks the next free port if 8321 is busy. |
# | **Large FITS** | Upload/mount data in cloud notebooks; preload with `files=`. |
# | **Two views** | Don't pass `open_browser=True` when embedding — you'll get a popup *and* an iframe. |
# | **Qt GUI** | `gui_qt()` needs a display server; prefer the web GUI in notebooks. |
# | **Stopping** | The server is a daemon thread; it stops when the kernel exits. `srv.stop()` is best-effort. |

# %%
# Detect environment (None outside IPython):
# mcf.detect_notebook_environment()
