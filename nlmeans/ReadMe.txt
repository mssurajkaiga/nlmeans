========================================
    Non-Local Means Denoiser Project
========================================

nlmeans.vcxproj
    This is the main project file for VC++ projects generated using an Application Wizard.
    It contains information about the version of Visual C++ that generated the file, and
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

nlmeans.vcxproj.filters
    This is the filters file for VC++ projects generated using an Application Wizard. 
    It contains information about the association between the files in your project 
    and the filters. This association is used in the IDE to show grouping of files with
    similar extensions under a specific node (for e.g. ".cpp" files are associated with the
    "Source Files" filter).

nlmeans.cpp
    This is the main application source file.

/////////////////////////////////////////////////////////////////////////////
Other standard files:

StdAfx.h, StdAfx.cpp
    These files are used to build a precompiled header (PCH) file
    named nlmeans.pch and a precompiled types file named StdAfx.obj.


/////////////////////////////////////////////////////////////////////////////

REMARKS:

General
=======
To maintain ease and flexibility in usage, the library is header-only with heavy templatization. Most part of the core library is optimized as much possible without breaking portability.

Non-Local Means Denoiser
========================
Implements patch-wise non-local means denoising as per Buades.et.al. The results are compared with those available on http://demo.ipol.im/demo/bcm_non_local_means_denoising/.
To view them, explore the results directory - running the matlab file will generate a basic comparison of the results.
You can either manually tweak the parameters or if you already have an estimate of the standard deviation (sigma) you can let the denoiser choose the rest of the parameters automatically.

Modified Non-Local Means Denoiser
=================================
Implements patch-wise dual-buffer non-local means denoising based on Rouselle.et.al. The input data can be generated in any monte-carlo render by splitting the samples in 2.
This was tested with a modified version of Mitsuba{Wenzel Jakob} which generated the necessary dual-buffer inputs.
An important point to note is that the scale of the vaariance cancellation parameter(alpha/sigma) mentioned in the paper is different from what is implemented here since it relies on the scale of the original
non-local means denoier configuration.
The parameters k and alpha here can radically affect the denoising results. Setting k too high results in dark patches whilst high values of alpha results in blocky artifacts. A favorable value would be in the range of 0.05f to 0.3f.
Unfortunately, there is no method yet to tweak the parameters automatically but perhaps the data from the variance estimate bitmaps could provide some knowledge on tweaking the settings {LOW PRIORITY}.