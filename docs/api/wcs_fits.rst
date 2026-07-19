WCS & FITS
==========

Headers, coordinates, cropping, and reprojection helpers.

.. currentmodule:: multicolorfits

Headers
-------

.. autosummary::
   :toctree: generated
   :nosignatures:

   McfHeader
   as_mcfheader
   save_header
   load_header
   force_header_2d
   force_header_3d
   force_header_floats
   make_simple_header
   squeeze_image
   annotate_provenance_header

Coordinates
-----------

.. autosummary::
   :toctree: generated
   :nosignatures:

   sex2dec
   dec2sex
   deg2dms
   dms2deg
   deg2hour
   hour2deg
   angular_distance
   sky_to_pixel
   pixel_to_sky
   header_coord_grids
   get_cdelts
   get_cd_matrix
   deg_per_pixel
   arcsec_per_pixel
   sterad_per_pixel
   beam_params_arcsec
   pixels_per_beam

Crop & reproject
----------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   crop_image
   crop_image_sky
   crop_cube
   crop_cube_sky
   reproject_image
   reproject_cube
   header_frame_name
   convert_header_frame
   reproject_to_frame
   reproject_to_galactic
   optimal_common_header

For long spectral cubes, :func:`~multicolorfits.reproject_cube` accepts
``print_progress=True`` (drives :func:`~multicolorfits.draw_progress_bar`).
See that helper's docstring for a ``tqdm`` alternative if you prefer wrapping
the per-plane loop yourself.
