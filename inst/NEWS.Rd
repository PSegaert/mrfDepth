\name{NEWS}
\title{News for package "mrfDepth"}
\encoding{UTF-8}

\section{Changes in version 1.0.4 (2017-06-13)}{
\itemize{
    \item FEATURE: DO or directional outlyingness now available
                   dirOutl, dprojdepth, dprojmedian
    \item FEATURE: functional measures now also have option for directional outlyingness. 
                   fOutl, mfd, mfdmedian
    \item BUGFIX: Call to warning updated in various functions, per new R call
 }
}

\section{Changes in version 1.0.3 (2017-04-01)}{
 \itemize{
    \item FEATURE: Added funtional depth functions: mfd, mfdmedian
    \item FEATURE: Added funtional outlyingness function: fOutl
    \item FEATURE: Added function distSpace
    \item FEATURE: Added the graphical tool FOM
    \item FEATURE: depthContour now also available for p = 1
    \item DOCUMENTATION: Various improvements to several help files
    \item DATA: Added the glass data set
    \item BUGFIX: Changed name of adjOutlyingness to adjOutl to avoid conflict 
                  with the robustbase package. 
    \item BUGFIX: Proper checking of scale for univariate bagdistance added 
    \item BUGFIX: Register native routines and disabled symbol search to fix 
                  warnings in R-devel
 } 
}

\section{Changes in version 1.0.2 (2017-02-01)}{
 \itemize{
    \item FEATURE: Added regression depth functions: rdepth, rdepthmedian and 
          cmltest.
    \item BUGFIX: Fix uninitialized use of variable in HalfSpaceDepthContour2D.f 
          to avoid warnings on fedora (as requested by CRAN).
 } 
}

\section{version 1.0.1 (2017-01-12)}{

 \itemize{
    \item Release on CRAN. 
 }
}