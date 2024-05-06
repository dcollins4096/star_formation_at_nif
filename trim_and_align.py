
from starter2 import *


#Coarse Cut from raw radiograph
#the import will extract the two frames from the raw radiograph.
#
# Produces the following images:
#     <frame>_full_images, the full tiff for <frame>
#     <frame>_noises, used to find the noise floor for extraction
#     <frame>_rough, the rough cut for each shot (by hand)
#     <frame>_trim, trim the image to the noise floor
import regions.regions_coarse as regions_coarse
trim_fname='data/trim.h5'
if not os.path.exists(trim_fname):
    trimmed = regions_coarse.do_coarse_trim()
    regions_coarse.save_trim(trimmed,trim_fname)
else:
    trimmed = regions_coarse.read_trim(trim_fname)

# Align the two shots to the point on the second fiducial.
# Produces:
#    <shot>_align1, the horizontal cuts through the fiducial.  Mostly a development tool
#    <frame>_trim_align, only the trimmed and aligned image
import regions.regions_align as regions_align
aligned_fname = 'data/aligned.h5'
regions_align.align_regions(trimmed, aligned_fname)


