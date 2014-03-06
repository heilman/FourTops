import os, sys, fileinput
import sendEmail
import commands

print 'Running complete tttt analysis'
print "############################################"
print "############################################"
print "############################################"
print "########           ELECTRON          #######"
print "########           CHANNEL           #######"
print "############################################"
print "############################################"
print "############################################"
os.system('rm log')
os.system('cp FourTop_EventSelection_El.cc FourTop_EventSelection_El.cc.back')

sendEmail.sendMail("Update from the fourtophawk","Starting full analysis run")

"""
########################################################################
#######               Run ScaleDown systematic                    ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
       sys.stdout.write(line.replace('int doScaleShift = 0;', 'int doScaleShift = 1;')) 
  
#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL > log')

print 'job completed successfully...'

#post-process
tail = "Scale down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)

########################################################################
#######               Run ScaleUp systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
	        sys.stdout.write(line.replace('int doScaleShift = 1;', 'int doScaleShift = 2;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log ')

tail = "Scale up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


#switch scale back to nominal
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
    sys.stdout.write(line.replace('int doScaleShift = 2;', 'int doScaleShift = 0;'))


########################################################################
#######               Run MatchingDown systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doMatchingShift = 0;', 'int doMatchingShift = 1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL > log')


tail = "Matching down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)
"""

########################################################################
#######               Run MatchingUp systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doMatchingShift = 0;', 'int doMatchingShift = 2;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "Matching up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)

for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                       sys.stdout.write(line.replace('int doMatchingShift = 2;', 'int doMatchingShift = 0;'))


########################################################################
#######               Run JESDown systematic                      ######
########################################################################

print "Running JES..."

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doJESShift = 0;', 'int doJESShift = 1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "JES down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)



########################################################################
#######               Run JESUp systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doJESShift = 1;', 'int doJESShift = 2;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "JES up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doJESShift = 2;', 'int doJESShift = 0;'))



########################################################################
#######               Run ttbbDown systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int dottbbShift = 0;', 'int dottbbShift = 1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "ttbb down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


########################################################################
#######               Run ttbbUp systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int dottbbShift = 1;', 'int dottbbShift = 2;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "ttbb up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int dottbbShift = 2;', 'int dottbbShift = 0;'))


########################################################################
#######               Run bTagDown systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int dobTagEffShift = 0;', 'int dobTagEffShift = -1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "btag down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)



########################################################################
#######               Run bTagUp systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int dobTagEffShift = -1;', 'int dobTagEffShift = 1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "btag up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int dobTagEffShift = 1;', 'int dobTagEffShift = 0;'))



########################################################################
#######               Run misTagDown systematic                   ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int domisTagEffShift = 0;', 'int domisTagEffShift = -1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "misTag down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)



########################################################################
#######               Run misTagUp systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int domisTagEffShift = -1;', 'int domisTagEffShift = 1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "misTag up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)

for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int domisTagEffShift = 1;', 'int domisTagEffShift = 0;'))




########################################################################
#######               Run LeptonSFDown systematic                 ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doLeptonSFShift = 0;', 'int doLeptonSFShift = 1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "leptonSF down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)




########################################################################
#######               Run LeptonSFUp systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doLeptonSFShift = 1;', 'int doLeptonSFShift = 2;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o EL')

#run
os.system('./EL >log')

#post-process
print 'job completed successfully...'

tail = "leptonSF up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doLeptonSFShift = 2;', 'int doLeptonSFShift = 0;'))




########################################################################
#######               Run PUDown systematic                       ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doPUShift = 0;', 'int doPUShift = 1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o El')

#run
os.system('./El >log')

#post-process
print 'job completed successfully...'

tail = "PU down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)



########################################################################
#######               Run PUUp systematic                         ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doPUShift = 1;', 'int doPUShift = 2;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o El')

#run
os.system('./El >log')

#post-process
print 'job completed successfully...'

tail = "PU up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doPUShift = 2;', 'int doPUShift = 0;'))




########################################################################
#######               Run JERDown systematic                      ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doJERShift = 0;', 'int doJERShift = 1;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o El')

#run
os.system('./El >log')

#post-process
print 'job completed successfully...'

tail = "JER down seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


########################################################################
#######               Run JERUp systematic                        ######
########################################################################

#edit alg to switch bool
for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doJERShift = 1;', 'int doJERShift = 2;'))


#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o El')

#run
os.system('./El >log')

#post-process
print 'job completed successfully...'

tail = "JER up seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


for i, line in enumerate(fileinput.input('FourTop_EventSelection_El.cc', inplace=1)):
                sys.stdout.write(line.replace('int doJERShift = 2;', 'int doJERShift = 0;'))
                



####Run nominal########

#compile
os.system('g++ -g -L ~/lib -I ../../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` FourTop_EventSelection_El.cc -o El')

#run
os.system('./El >log')

#post-process
print 'job completed successfully...'

tail = "Nominal seems to have completed successfully      " + str(commands.getstatusoutput('tail log'))

sendEmail.sendMail("Update from the fourtophawk",tail)


print "############################################"
print "############################################"
print "############################################"
print "########           END               #######"
print "########           JOB               #######"
print "############################################"
print "############################################"
print "############################################"
