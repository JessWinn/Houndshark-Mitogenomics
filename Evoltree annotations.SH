## Evolview Tree Design ##

##style 1
##!grouplabel	color=pink,fontsize=14,fontcolor=black,linewidth=2
##style 2
##!grouplabel	style=2,color=pink,show=1,marginPCT=0.05,fontsize=14,fontcolor=white,fontitalic=0,textalign=middle,textorientation=horizontal,linewidth=2
##style 3
##!grouplabel	style=3,color=pink,show=1,marginPCT=0.05,fontsize=14,fontcolor=white,fontitalic=0,textalign=middle,textorientation=horizontal,linewidth=2
##style 4
##!grouplabel	style=4,color=pink,show=1,marginPCT=0.05,fontsize=14,fontcolor=white,fontitalic=0,textalign=middle,textorientation=horizontal,linewidth=2
##style 5
##!grouplabel	style=5,color=pink,show=1,marginPCT=0.05,fontsize=14,fontcolor=white,fontitalic=0,textalign=middle,textorientation=horizontal,linewidth=2
!op	0.8

##style 1	leaf_name	text=mammal,color=darkgreen,textorientation=vertical,linewidth=4,fontsize=16,linestyle=dashed
##style 2-5	leaf_name	bkcolor=#BE4144,text=mammal,textorientation=vertical,linewidth=4,fontsize=16
!grouplabel	style=2
!op	0.8
Carcharhinus_brachyurus_NC_057525,Carcharhinus_leucas_NC_023522	bkcolor=#FCDABE,text=Carcharhinidae
Sphyrna_lewini_NC_022679,Sphyrna_zygaena_NC_025778	bkcolor=#B7F1A5,text=Sphyrnidae
Hemigaleus_microstoma_NC_029400,Hemipristis_elongata_NC_032065	bkcolor=#FF9AA2,text=Hemigaleidae
Hemitriakis_japanica_NC_026774,Galeorhinus_galeus_ON652874	bkcolor=#B2EBF9,text=Triakidae
Proscyllium_habereri_NC_030216	bkcolor=#FFFEC4,text=Proscyllidae
Pseudotriakis_microdon_NC_022735	bkcolor=#F9D2EF,text=Pseudotriakidae
Cephaloscyllium_umbratile_NC_029399,Poroderma_pantherinum_NC_043830	bkcolor=#C6B7F1,text=Scyliorhinidae
Carcharodon_carcharias_NC_022415	bkcolor=lightgrey,text=Outgroup

# Mustelus genus reproduction mode #

!grouplabel	style=1
!op	0.8
Mustelus_asterias_ON652873,Triakis_megalopterus_ON075075	text=aplacental-spotted,color=#0033CC,textorientation=horizontal,linewidth=4,fontsize=14,linestyle=dashed
Mustelus_griseus_NC_023527,Mustelus_mustelus_NC_039629	text=placental-non-spotted,color=#3399FF,textorientation=horizontal,linewidth=4,fontsize=14 

# Bootstrap style #

!bootstrapValueStyle	show=2,style=numeric,color=darkred,place=2,size=12

# Marking our sequences #

##single decoration shape
##double decoration shapes seperated by tab space
##leaf_name	shape,fill_color
##leaf_name	shape,fill_color:border_color
Mustelus_asterias_ON652873	star,red	
Mustelus_palumbes_ON075076	star,red
Triakis_megalopterus_ON075075	star,red
Mustelus_mosis_ON075077	star,red
Galeorhinus_galeus_ON652874	star,red
