\name{NEWS}
\title{News for package "mrfDepth"}
\encoding{UTF-8}

\section{Changes in version 1.0.10 (2018-10-12)}{
\itemize{
    \item Bugfixes in distSpace method
    \item Update of references in documentation
 }
}

\section{Changes in version 1.0.9 (2018-09-26)}{
\itemize{
    \item Change of maintainer.
    \item Change of cutoff in outlyingness, adjOutl and dirOutl
    \item Change of compiler flags (removed C flags)
 }
}


\section{Changes in version 1.0.8 (2018-07-08)}{
\itemize{
    \item BUGFIX: Componentwise option of Directional Outlylingness fixed
    \item Feature: Weighted sd implemented in FOM
    \item MISC: Cleanup of C++ code for Directional Outlyingness
 }
}




\section{Changes in version 1.0.7 (2018-07-01)}{
\itemize{
    \item Feature: Added functional data examples of MFOD / MFC paper: 
                   characterA,
                   characterI,
                   mri,
                   plane,
                   tables,
                   wine
                   
    \item MISC: Documentation of random number generator 
                for normal distribution used by FORTRAN code.
 }
}


\section{Changes in version 1.0.6 (2017-11-20)}{
\itemize{
    \item BUGFIX: Changed internal call
                  compiler:::cmpfun
                  to 
                  compiler::cmpfun
                  as requested by CRAN.
 }
}


\section{Changes in version 1.0.5 (2017-10-15)}{
\itemize{
    \item BUGFIX: Bugfix for Solaris as requested by CRAN.
                  -> Fixing "dirOutl.h:85:55: error""
    \item BUGFIX: Bugfix for Windows and Solaris as requested by CRAN.
                  -> Makevars.win file
    \item BUGFIX: Corrections of typos in manuals 
 }
}


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