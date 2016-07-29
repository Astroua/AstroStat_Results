
'''
Script for reproducing the figures in Appendix A

Comparisons are between Designs 2 & 24.
'''

from turbustat.data_reduction import Mask_and_Moments
from turbustat. statistics import Wavelet_Distance, MVC_Distance, \
    PSpec_Distance, BiSpectrum_Distance, GenusDistance, \
    DeltaVariance_Distance, VCA_Distance, VCS_Distance, Tsallis_Distance, \
    StatMoments_Distance, PCA_Distance, SCF_Distance, Cramer_Distance, \
    DendroDistance, PDF_Distance
import os
import matplotlib.pyplot as p
import seaborn as sns
sns.set_style("white")
import astropy.units as u

p.ioff()


path_to_data = "/media/eric/Data_3/SimSuite8/"
moments_path = os.path.join(path_to_data, "moments/")


des2 = "lustrehomeerosSimSuite8Design2_flatrho_0029_00_radmc.fits"
dataset1 = Mask_and_Moments.from_fits(os.path.join(path_to_data, des2),
                                      moments_path=moments_path).to_dict()
# des19 = "lustrehomeerosSimSuite8Design19_flatrho_0030_00_radmc.fits"
des19 = "lustrehomeerosSimSuite7Fiducial1_flatrho_0029_00_radmc.fits"
dataset2 = Mask_and_Moments.from_fits(os.path.join(path_to_data, des19),
                                      moments_path=moments_path).to_dict()

label1 = "Design 2"
# label2 = "Fiducial 1"
label2 = "Design 19"

values = {}

# vcs_distance = VCS_Distance(dataset1["cube"],
#                             dataset2["cube"], breaks=-0.5)
# vcs_distance.distance_metric(verbose=True, label1=label1, label2=label2)

# Wavelet Transform

wavelet_distance = \
    Wavelet_Distance(dataset1["integrated_intensity"],
                     dataset2["integrated_intensity"])
wavelet_distance.distance_metric(verbose=True, label1=label1, label2=label2)

print "Wavelet Distance: %s" % (wavelet_distance.distance)
# values["Wavelet_Slopes"] = [wavelet_distance.]

# p.savefig()

# MVC

mvc_distance = MVC_Distance(dataset1, dataset2)
mvc_distance.distance_metric(verbose=True, label1=label1, label2=label2)

print "MVC Distance: %s" % (mvc_distance.distance)

# p.savefig()

# Spatial Power Spectrum/ Bispectrum

pspec_distance = \
    PSpec_Distance(dataset1["integrated_intensity"],
                   dataset2["integrated_intensity"],
                   weights1=dataset1["integrated_intensity_error"][0]**2.,
                   weights2=dataset2["integrated_intensity_error"][0]**2.)
pspec_distance.distance_metric(verbose=True, label1=label1, label2=label2)

print "Spatial Power Spectrum Distance: %s" % (pspec_distance.distance)

# p.savefig()

bispec_distance = \
    BiSpectrum_Distance(dataset1["integrated_intensity"],
                        dataset2["integrated_intensity"])
bispec_distance.distance_metric(verbose=True, label1=label1, label2=label2)

print "Bispectrum Distance: %s" % (bispec_distance.distance)

# p.savefig()

# Genus

genus_distance = \
    GenusDistance(dataset1["integrated_intensity"][0],
                  dataset2["integrated_intensity"][0])
genus_distance.distance_metric(verbose=True, label1=label1, label2=label2)

print "Genus Distance: %s" % (genus_distance.distance)

# p.savefig()

# Delta-Variance

delvar_distance = \
    DeltaVariance_Distance(dataset1["integrated_intensity"],
                           dataset2["integrated_intensity"],
                           weights1=dataset1["integrated_intensity_error"][0],
                           weights2=dataset2["integrated_intensity_error"][0])
delvar_distance.distance_metric(verbose=True, label1=label1, label2=label2)

print "Delta-Variance Distance: %s" % (delvar_distance.distance)

# p.savefig()

# VCA/VCS

vcs_distance = VCS_Distance(dataset1["cube"],
                            dataset2["cube"], breaks=-0.5)
vcs_distance.distance_metric(verbose=True, label1=label1, label2=label2)

print "VCS Distance: %s" % (vcs_distance.distance)

# p.savefig()

vca_distance = VCA_Distance(dataset1["cube"],
                            dataset2["cube"])
vca_distance.distance_metric(verbose=True, label1=label1, label2=label2)

print "VCA Distance: %s" % (vca_distance.distance)

# p.savefig()

# Tsallis

tsallis_distance = Tsallis_Distance(dataset1["integrated_intensity"][0],
                                    dataset2["integrated_intensity"][0])
tsallis_distance.distance_metric(verbose=True)

# p.savefig()

print "Tsallis Distance: %s" % (tsallis_distance.distance)

# High-order stats

moment_distance = StatMoments_Distance(dataset1["integrated_intensity"][0],
                                       dataset2["integrated_intensity"][0], 5)
moment_distance.distance_metric(verbose=True, label1=label1,
                                label2=label2)

print "Kurtosis Distance: %s" % (moment_distance.kurtosis_distance)

print "Skewness Distance: %s" % (moment_distance.skewness_distance)

# p.savefig()

# PCA

pca_distance = PCA_Distance(dataset1["cube"][0],
                            dataset2["cube"][0],
                            normalize=False, mean_sub=True)
pca_distance.distance_metric(verbose=True, label1=label1,
                             label2=label2)

print "PCA Distance: %s" % (pca_distance.distance)

# p.savefig()

# SCF

scf_distance = SCF_Distance(dataset1["cube"],
                            dataset2["cube"], size=21)
scf_distance.distance_metric(verbose=True, label1=label1,
                             label2=label2)

print "SCF Distance: %s" % (scf_distance.distance)

# p.savefig()

# Dendrogram Stats

dendro_distance = \
    DendroDistance(dataset1["cube"][0],
                   dataset2["cube"][0])
dendro_distance.distance_metric(verbose=True, label1=label1,
                                label2=label2)

print dendro_distance.num_distance
print dendro_distance.histogram_distance

# p.savefig()

# PDF

pdf_distance = \
    PDF_Distance(dataset1["integrated_intensity"][0],
                 dataset2["integrated_intensity"][0],
                 min_val1=0.0,
                 min_val2=0.0,
                 weights1=dataset1["integrated_intensity_error"][0] ** -2.,
                 weights2=dataset2["integrated_intensity_error"][0] ** -2.)

pdf_distance.distance_metric(verbose=True, show_data=False,
                             label1=label1, label2=label2)

print("Hellinger Distance: " + str(pdf_distance.hellinger_distance))
print("KS Distance: "+str(pdf_distance.ks_distance),
      "KS p-value"+str(pdf_distance.ks_pval))

# p.savefig()
