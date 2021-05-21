import os

def xfm_subcortToHcp(ants_path, freesurfer_path, wb_path, recon_all_folder, segmentation_folder, subcortex_label, mni2mm, empty_cifti_cortex_left, empty_cifti_cortex_right, template_cifti, TR, workdir, outputdir):
    
    '''
    This script maps volumetric subThalamic segmentations to cifti volume to 
    make a mask. Some of the inputs can be find in the utilities folder of this 
    repo (in the segmentTocifti_files.zip archive).
    
    Inputs:         
        recon_all_folder = Unzipped folder containing the main recon-all results 
        segmentation_folder = Unzipped register_thalamic_nuclei gear output folder 
        subcortex_label = label for creating subcortex
        mni2mm = MNI template 2mm
        empty_ciftify_cortex_left = left cortex with all zeros
        empty_ciftify_cortex_right = right cortex with all zeros
        template_cifti = template cifti dtseries
        TR = TR in seconds
        workdir = workdir where intermediate files will be saved 
        outputdir = outputdir where the final map will be saved 
    '''
    
    # Create folders if don't exist
    if not os.path.exists(workdir):
        os.system('mkdir %s' % workdir)
    if not os.path.exists(outputdir):
        os.system('mkdir %s' % outputdir)
        
    # Register orig to MNI
    orig_file_mgz = os.path.join(recon_all_folder, 'mri', 'orig.mgz')
    orig_file_nii = os.path.join(workdir, 'orig.nii.gz')
    os.system('%s %s %s' % (os.path.join(freesurfer_path, 'mri_convert'), orig_file_mgz, orig_file_nii))
    
    antsReg = ('%s --verbose 1 --dimensionality 3 --float 0 \
--collapse-output-transforms 1 --output [ %s/orig2MNI,%s/orig2MNIWarped.nii.gz,%s/orig2MNIInverseWarped.nii.gz ] \
--interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] \
--initial-moving-transform [%s,%s,1] \
--transform Rigid[ 0.1 ] --metric MI[%s,%s,1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
--transform Affine[ 0.1 ] --metric MI[%s,%s,1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
--transform SyN[0.1,3,0] --metric CC[%s,%s,1,4] \
--convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox' % (os.path.join(ants_path, 'antsRegistration'), workdir, workdir, workdir, mni2mm, orig_file_nii,
                                                                                                mni2mm, orig_file_nii, mni2mm, orig_file_nii,
                                                                                                mni2mm, orig_file_nii))
    os.system(antsReg)
    generic_affine = os.path.join(workdir, 'orig2MNI0GenericAffine.mat')
    warp = os.path.join(workdir, 'orig2MNI1Warp.nii.gz')
    
    # Move segmentations to MNI 
    segmentation_mgz = os.path.join(segmentation_folder, 'ThalamicNuclei.v12.T1.FSvoxelSpace.mgz')
    segmentation_nifti = os.path.join(workdir, 'ThalamicNuclei.v12.T1.FSvoxelSpace.nii.gz')
    os.system('%s %s %s' % (os.path.join(freesurfer_path, 'mri_convert'), segmentation_mgz, segmentation_nifti))
    
    segmentations_in_mni = os.path.join(workdir, 'MNI_ThalamicNuclei.v12.T1.FSvoxelSpace.nii.gz')
    os.system('%s -d 3 -i %s -r %s -t %s -t %s -n NearestNeighbor -o %s' % (os.path.join(ants_path, 'antsApplyTransforms'), segmentation_nifti, mni2mm, 
                                                                            warp, generic_affine, segmentations_in_mni))
    
        
    # Make the cifti label file
    os.system('%s -cifti-create-dense-from-template %s %s -series %s 0 -volume-all %s -metric CORTEX_LEFT %s -metric CORTEX_RIGHT %s' % (os.path.join(wb_path, 'wb_command'), template_cifti, os.path.join(outputdir, 'subcortex_cifti.dtseries.nii'),
                                                                                                                                         TR, segmentations_in_mni, empty_cifti_cortex_left, empty_cifti_cortex_right))
    
    return segmentations_in_mni
