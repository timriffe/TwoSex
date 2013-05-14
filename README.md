"The two-sex problem in populations structured by remaining years of life" 
----------------------------------------------------------------
Welcome to my dissertation code dump.
For the time being, this is the only route to reproducing my results and I apologize if it's a labyrinth!

This includes ALL R scripts and tex files
----------------------------------------------------------------
It does NOT include three folders, which are not tracked, and which therefore do not appear here, but that 
are also integral parts of the workflow:

/Data/             (files are nearly all in .Rdata format, so not good to track with git. Download separately as stated above)

/DiagnosticPlots/  (these are not in the dissertation. These figures were early looks at the data)

/latex/Figures/    (folder to hold pdfs for inclusion in the main latex document. You can reproduce 
                    its contents with the R scripts and data)
                    
All figures can be reproduced by using the R script files if you have the data, 
which you need to download separately from Google Drive under the following link:

https://drive.google.com/folderview?id=0B2b7NmVR77Q1bDdSNkVaVGxJZVE&usp=sharing

This will give you a folder called 'Data', which should be placed in /DISS/ 
at the same level as /R/  (and friends)

Once you have the data, R scripts will (hopefully) be able to grab it and use it to make
Figures, or to help you test functions.
All R functions are undocumented, though some are lightly annotated (but not thoroughly at this time). 
I'll be trickle-grooming this code, so either make a request or be patient. 
In general, you'll need to tinker to make sense of things, sorry.

Each R script begins by setting the working directory with, with setwd(), 
so you need to change that appropriately

*DO NOT* mindlessly run an entire script without looking at it first, as it may contain a simulation
that will dominate your computer resources.

thesisMain.tex is the mother tex document, pulling in chunks from /latex/ 
The structure of \thesisMain.tex is what gives the dissertation its part/chapter/section structure
If you want to find how functions were implemented, the best way (for now) is to either 1) intuit via
script names, or 2) to find a figure in the dissertation that used this code, then find the .tex file that
makes that section. In the latex chunk where the figure is included you'll (usually) see the name of the script
that produced it.

Defense is on June 26th, 2013 UAB, Barcelona. 
---------------------------------------------
wish me luck

