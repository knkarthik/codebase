#patients is a dataframe with all the information of course.
#Make a 3-way frequency table
arr = table(patients$Patient, patients$Familiarl.sporadic, patients$ACMG.classification)

#COnvert it into array
arr = (unclass(arr))

#Add dimnames
dimnames(arr) = list(patients = unique(patients$Patient), History=c("Familial", "sporadic"), Classification=c("Benign ", "Likely benign","Likely pathogenic","Pathogenic ","Uncertain significance"))

# #str(arr) will look like:
# int [1:10, 1:2, 1:5] 0 0 0 0 0 0 0 0 0 1 ...
# - attr(*, "dimnames")=List of 3
# ..$ patients      : chr [1:10] "7" "9" "46" "60" ...
# ..$ History       : chr [1:2] "Familial" "sporadic"
# ..$ Classification: chr [1:5] "Benign " "Likely benign" "Likely pathogenic" "Pathogenic " ...

#Make a frequency table
ftable(arr, row.vars=c("History","Classification"), col.vars="patients")
#should look something like:
#                                  patients 7 9 46 60 73 76 84 92 127 128
# History   Classification                                               
# Familial  Benign                          0 0  0  0  0  0  0  0   0   1
#           Likely benign                   2 0  0  0  0  0  0  0   0   1
#           Likely pathogenic               0 0  0  0  0  2  0  0   0   0
#           Pathogenic                      0 0  0  0  0  0  0  0   0   0
#           Uncertain significance          4 0  0  0  0  3  0  0   0   4
# sporadic  Benign                          0 1  0  0  0  0  0  0   0   0
#           Likely benign                   0 3  2  1  1  0  2  1   1   0
#           Likely pathogenic               0 0  0  0  0  0  0  0   0   0
#           Pathogenic                      0 0  0  0  0  0  0  1   0   0
#           Uncertain significance          0 2  0  4  1  0  1  1   2   0