# TSUFLIND_pd

All the parameter necessary to run TSUFLIND are in the file parameter_P14abc.txt. The filename does not make much sense, but please don't change it. If you open the file in a text editor, you will see that there are 59 lines in that file, and every line not containing # is important. However, most of the lines do not need to be changed and most of the parameters are self explanatory. The first important line is 37: water_run_up, which is exactly what you have. Line 39 is also important, which contains the slope given as 1/. Lines 37 and 39 are both used in the Soulsby model. If you do not know these values, please estimate them. While it is good to know these parameters in great detail, I will provide more information on that topic in a later blog post. Line 43 denotes the number of samples you have and line 47 gives the file name (sample_P14abc.csv) containing the respective grain-size distributions of the samples. The last important line is 53. The named file (test.csv) contains the location along the slope and the total thicknesses of the deposits. 

Running TSUFLIND and Final Remarks
When you are done with your modifications of parameter_P14abc.txt and test.csv you will need to run main.py. From the command line you can type python main.py. If you see something like the following..you are good:
bash prompt$ python main.py 
# Outputting to File[data_process.csv]:
[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
# Outputting to File[sample00.csv]:
[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
# Outputting to File[sample01.csv]:
[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
and then:
====================================running==================================== 
...everything is fine. Just give it some time to do the magic.
Let me make some final remarks: TSUFLIND is a research code at is present stage and this not a full manual. However, it gives a software to "play" around and learn about sensitivities. Playing around with models is extremely important in order to learn what they can and cannot do. If you did not change any the files ending with py and the model does not run you did not screw up the model...just the input files. I hope to write another blog post about TSUFLIND soon that will serve more as a tutorial. If you have any questions, trouble, or want to use TSUFLIND for a specific example that you cannot get to run, please to not hesitate to contact me at weiszr@vt.edu
