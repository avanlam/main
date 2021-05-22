function rank_grp = group_ranking ( rank_indvl, numGrp )
numDim = length (rank_indvl);
grpSize = round ( numDim / numGrp );
grpNum = 1; temp = 0;
for rankNum = 1 : numDim
    if temp == grpSize
        temp = 0;
        grpNum = grpNum + 1;
    end
    temp = temp + 1;
    rank_grp ( 1,  rank_indvl == rankNum  ) = grpNum;
end
end