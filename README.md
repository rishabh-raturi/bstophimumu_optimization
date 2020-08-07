*  **STEP 1**

Get the value of effective sigma from rare mode MC fit (no cuts applied). Sigma is used to define the signal and side-band region.

**Use script-> dcb.C**

This script uses double crystal ball with twp tail parameters- alpha1 and alpha2. Give positive range to one and negative to another.  

*  **STEP 2**

Plot the signal side-band overlaying plots to find the discriminating variables used as inputs. 

For the current analysis the input variables used are-
1.  B(s) vtx. CL
2.  B(s) cosalpha 2D
3.  B(s) L_{xy}
4.  TMath::Max(kaons pt)
5.  TMath::Max(muons pt)
6.  TMath::Max(kaons trk. MinIP)


*  **STEP 3**

Create reduced ntuples for BDT training and testing. 

**Use scripts->**
*  selectSIGNAL.cc
*  selectBKG.cc

Use extended side-band for training to have more statistics in data, which will improve the performance of the BDT classifier


*  **STEP 4**


Train BDT

**USe script-> TMVAClassification.C**

Output-> in the form of GUI, save the important plots like overtraining plot, input variables, co-relation matrix, ROC curve.


*  **STEP 5**

Best B(s) selection based on the best BDT score.


*  **STEP 6**

After step 5 you have the ntuples which has the BDT classifier stored in it, plot the significance cutve to get the optimized BDT cut value. 

**Use script-> BDT_cut_S_over_SplusB.cc**


*  **STEP 7**

Create reduced ntuples for fitting. Cuts need to be applied are mentoned in the scripts.

Use the scripts stored in the directory-> **reduced_ntuples**




*  **Step 8**

Perform MC fits for all the modes and then use the parameters to fit data. 

