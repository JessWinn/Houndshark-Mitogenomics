## STEP 8: annotating species trees in Evolview v.3 webserver (https://www.evolgenius.info/evolview/)

1. Open the treefiles for Maximum Likelihood PS5 and ASTRAL (elasmo-mitophy-15G-ASTRAL.tre) in FigTree v.1.4.4. and add bootstrap annotations. 

Save the file with annotations as a NEX file.

2. Create a new project and open the NEX file in Evolview.

3. Add the following annotations.

# Orders

!grouplabel	style=1
!op	0.8
NC_023522.1_Carcharhinus_leucas,NC_043830.1_Poroderma_pantherinum	text=Carcharhiniformes,color=#27282A,textorientation=vertical,linewidth=4,fontsize=12
NC_024269.1_Lamna_ditropis,NC_022822.1_Alopias_pelagicus	text=Lamniformes,color=#626469,textorientation=vertical,linewidth=4,fontsize=11
KC633221.1_Rhincodon_typus,NC_022148.1_Orectolobus_japonicus	text=Orectolobiformes,color=#A4A5A9,textorientation=vertical,linewidth=4,fontsize=11

## Families ##

!grouplabel	style=1
!op	0.8
NC_023522.1_Carcharhinus_leucas,NC_029843.1_Loxodon_macrorhinus	text=Carcharhinidae,color=#008EED,textorientation=horizontal,linewidth=4,fontsize=14 
NC_022193.1_Galeocerdo_cuvier	text=Carcharhinidae, color=#008EED,textorientation=horizontal,linewidth=4,fontsize=14 
NC_022679.1_Sphyrna_lewini,NC_031812.1_Eusphyra_blochii	text=Sphyrnidae, color=#FF6821, textorientation=horizontal,linewidth=4,fontsize=14 
NC_029400.1_Hemigaleus_microstoma,NC_032065.1_Hemipristis_elongata	text=Hemigaleidae, color=green, textorientation=horizontal,linewidth=4,fontsize=14 
KJ617039.1_Hemitriakis_japanica,ON652874_Galeorhinus_galeus	text=Triakidae,color=red, textorientation=horizontal,linewidth=4,fontsize=14 
NC_022735.1_Pseudotriakis_microdon	text=Pseudotriakidae, color=purple,textorientation=horizontal,linewidth=4,fontsize=14
NC_030216.1_Proscyllium_habereri	text=Proscyllidae,color=#EAC400,textorientation=horizontal,linewidth=4,fontsize=14 
NC_031811.1_Halaelurus_buergeri,NC_049881.1_Galeus_melastomus	text=Pentanchidae,color=turquoise,textorientation=horizontal,linewidth=4,fontsize=14
NC_056784.1_Parmaturus_melanobranchus,NC_043830.1_Poroderma_pantherinum	text=Scyliorhinidae,color=#F4A460,textorientation=horizontal,linewidth=4,fontsize=14

# Mustelus reproduction modes

##default style is 1
##style 2 is gradient based
!leaf_annotation_style	2

!groups_class	id=T1,text=,textalign=middle,shape=rect,bkcolor=lightpink,color=black,fontsize=12
!groups_class	id=T2,text=,textalign=middle,shape=rect,bkcolor=powderblue,color=black,fontsize=12

ON652873_Mustelus_asterias,ON075075_Triakis_megalopterus	T1
NC_023527.1_Mustelus_griseus,NC_039629.1_Mustelus_mustelus	T2

# Our sequences

ON652873_Mustelus_asterias	red
ON075076_Mustelus_palumbes	red
ON075075_Triakis_megalopterus	red
ON075077_Mustelus_mosis	red
ON652874_Galeorhinus_galeus	red	

# Bootstraps

!bootstrapValueStyle	show=1,style=circle,place=3,strokewidth=3
0,94	color=white:blue
95,100	color=none
