
from starter2 import *
import regions.regions_coarse as regions_coarse
import regions.regions_align as regions_align
reload(regions_align)
reload(regions_coarse)


#Coarse Cut from raw radiograph
#the import will extract the two frames from the raw radiograph.
#
# Produces the following images:
#     <frame>_full_images, the full tiff for <frame>
#     <frame>_noises, used to find the noise floor for extraction
#     <frame>_rough, the rough cut for each shot (by hand)
#     <frame>_trim, trim the image to the noise floor
trim_fname='data/trim.h5'
if not os.path.exists(trim_fname):
    trimmed = regions_coarse.do_coarse_trim()
    regions_coarse.save_trim(trimmed,trim_fname)
else:
    trimmed = regions_coarse.read(trim_fname)

# Align the two shots to the point on the second fiducial.
# Produces:
#    <shot>_align1, the horizontal cuts through the fiducial.  Mostly a development tool
#    <frame>_trim_align, only the trimmed and aligned image
aligned_fname = 'data/aligned.h5'
if not os.path.exists(aligned_fname) or True:
    #regions_align.align_regions(trimmed, aligned_fname)
    #regions_align.align_regions2(trimmed, align_fname=aligned_fname, shot_list=['s120'])
    regions_align.align_regions2(trimmed, align_fname=aligned_fname, shot_list=['r0','r60','r120','s90','s120'])


