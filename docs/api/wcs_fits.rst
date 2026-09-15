WCS & FITS
==========

Headers, coordinates, cropping, and reprojection helpers.

* Guide: :doc:`/guide/wcs_and_alignment`, :doc:`/guide/preparing_images`
* Job snippets: :doc:`/capabilities/index`
* Stack align / prep also listed under :doc:`pipeline`
  (``align_stack``, ``prep_layers``)

.. currentmodule:: multicolorfits

Headers
-------

``describe_header(..., verbose=False)`` returns the summary dict without
printing. Stretch suggestions live on :func:`~multicolorfits.describe_image`
/ :func:`~multicolorfits.describe_images` (:doc:`/guide/intensity_scaling`).
Start everyday scripts with :func:`~multicolorfits.tidy_header` /
:func:`~multicolorfits.read_fits`.

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
   tidy_header
   describe_header
   wcs_is_flipped
   east_increases_right
   squeeze_image
   annotate_provenance_header
   read_fits
   write_fits

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
   convolution_beam
   match_beam
   match_beam_to_header
   convolve2Dgaus
   convolve2Dgaus_matchhdr

Crop & reproject
----------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   crop_image
   crop_image_sky
   crop_cube
   crop_cube_sky
   crop_to_overlap
   blank_missing
   reproject_image
   reproject_cube
   header_frame_name
   convert_header_frame
   reproject_to_frame
   reproject_to_galactic
   make_rotated_header
   make_north_up_header
   reproject_to_rotation
   reproject_north_up
   optimal_common_header

For long spectral cubes, :func:`~multicolorfits.reproject_cube` accepts
``print_progress=True`` (drives :func:`~multicolorfits.draw_progress_bar`).
See that helper's docstring for a ``tqdm`` alternative if you prefer wrapping
the per-plane loop yourself.
