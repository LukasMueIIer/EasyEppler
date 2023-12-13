import numpy as onp

class AirFoil:
    def __init__(self,_FileName:str,_Initials:str,_Number:int,_N:int = 60,_ZeroN = 31):
        self.N = _N #amount of points the foil has
        self.ZeroN = _ZeroN #to compensate for the mapping not being an exact circle -> the N value that closly coincides with LE 
        #!! must be before LE when coming from the TOP starting from TE
        self.FileName = _FileName #File in which instructions will be written
        self.Initials = _Initials #Leading Letters of foil Name
        self.Number = _Number     #Airfoil Number
        self.UseRamp = 0          #Checks if ramp parameter have been set
        self.dc_Front_Up = 0
        self.dc_Back_Up = 0
        self.dc_Front_Low = 0
        self.dc_Back_Low = 0
        
        self.index_RE = 0
        self.REs = onp.zeros(5) #RE numbers for visous calculations
        self.index_transMode = onp.zeros(5)
        self.transTop = onp.zeros(5)
        self.transBottom = onp.zeros(5)
    
    def CtoNUpper(self,c:float) -> float: #converts a chord position on the upper side to N
        if(c == 0):
            print("WARNING using c=0 unless at LE will CRASH eppler!!!!!")
        elif(c < 0):
            print("WARNING c can't be smaller than 0")
        elif(c > 1):
            print("WARNING c can't be bigger than 1")
        m = onp.arccos(2*c - 1) * (self.N / 2) / onp.pi
        m = m/(0.5 * self.N) * self.ZeroN #correction for not exact circle
        return m

    def CtoNLower(self,c:float) -> float: #converts a chord position on the lower side to N
        if(c == 0):
            print("WARNING using c=0 unless at LE will CRASH eppler!!!!!")
        elif(c < 0):
            print("WARNING c can't be smaller than 0")
        elif(c > 1):
            print("WARNING c can't be bigger than 1")
        m = onp.arccos(1 - 2*c) * (self.N / 2) / onp.pi
        m = m/(0.5 * self.N) * (self.N - self.ZeroN) #correction for not exact circle
        m = m + self.ZeroN
        return m

    def CommandSyntaxHelper(self,Command:str,NUPA,NUPE,NUPI,NUPU):
        _string = Command + str(NUPA) + str(NUPE) + str(NUPI) + str(NUPU)
        if(len(_string) != 10):
            print("WARNING Command String size is wrong!!!")
        self.file.write(_string + " ")

    def WriteUpperC_Alpha(self,c,a) -> int: #writes an array of cpositions and their respektive alpha where they should be constant to TRA1 Line
        #convert c's to nues and check that they increase monotonically
        nue = c
        nuebefore = -0.001
        for i in range(0,len(c)):
            nue[i] = self.CtoNUpper(c[i])
            if(nue[i] <= nuebefore):
                print("WARNING Nue must be increasing")
            nuebefore = nue[i]
        
        if(len(c) != len(a)):
            print("WARNING arrays are not same size")
            return 1

        File = self.File

        File.write("TRA1" + " " + " " + str(self.Number).zfill(4)) #begin TraLine
        
        for i in range(0,len(c)-1):
            if(i > 5):
                if((i % 10) == 0):
                    File.write("\n")
                    File.write("TRA1" + " " + " " + str(self.Number).zfill(4)) #begin new TraLine we dont want more than 10 entries
            File.write(" " + str(round(nue[i],1)) + " " + str(round(a[i],1))) 

        #check for last digit
        if(self.ZeroN - nue[-1] < 0.3):
            File.write(" " + str(round(0,1)) + " " + str(round(a[-1],1)))
        else:
            File.write(" " + str(round(nue[-1],1)) + " " + str(round(a[-1],1)))
            File.write(" " + str(round(0,1)) + " " + str(round(a[-1],1)))  #we need a zero entrie
            print("WARNING Zero entrie was automatticaly added")
        File.write("\n")
        return 0
    
    def WriteLowerC_Alpha(self,c,a,midSnap = 1) -> int: #writes an array of cpositions and their respektive alpha where they should be constant to TRA1 Line
        #Mid snap moves points to the closest .5 number accodring to Eppler this yeilds smother Airfoils
        #convert c's to nues and check that they increase monotonically
        nue = c
        nuebefore = -0.001
        for i in range(0,len(c)):
            nue[i] = self.CtoNLower(c[i])
            if(nue[i] <= nuebefore):
                print("WARNING Nue must be increasing")
            nuebefore = nue[i]
        
        if(len(c) != len(a)):
            print("WARNING arrays are not same size")
            return 1

        File = self.File

        File.write("TRA1" + " " + " " + str(self.Number).zfill(4)) #begin TraLine
        
        for i in range(0,len(c)-1):
            if(i > 5):
                if((i % 10) == 0):
                    File.write("\n")
                    File.write("TRA1" + " " + " " + str(self.Number).zfill(4)) #begin new TraLine we dont want more than 10 entries
            File.write(" " + str(round(nue[i],1)) + " " + str(round(a[i],1))) 

        #check for last digit
        if(self.N - nue[-1] < 0.51):
            File.write(" " + str(round(self.N,1)) + " " + str(round(a[-1],1)))
        else:
            File.write(" " + str(round(nue[-1],1)) + " " + str(round(a[-1],1)))
            File.write(" " + str(round(self.N,1)) + " " + str(round(a[-1],1)))  #we need a Final Entry
        File.write("\n")
        return 0

    def BeginnFile(self) -> int: #Write all the stuff to file
        self.File = open(self.FileName, "w")
        File = self.File
        File.write("REMO1       *8" + self.Initials + "\n")
        File.write("REMO1      *P@1A@2IRFOIL DESIGN 2020/2021" + "\n")
        #Park ich mal hier print(str(1).zfill(2)) Für später Airfoil Number
        return 0

    def MPRUpper(self,cMPR,cClosure,par1,par2,mode=2):
        if(cClosure < cMPR):
            print("WARNING cMPR shouldn't be bigger than cClosure")
        self.MPRUpNue = self.CtoNUpper(cMPR)
        self.cMPRup = cMPR
        self.CloseUpNue = self.CtoNUpper(cClosure)
        self.MPRModeUper = mode
        self.MPRF1Uper = par1
        self.MPRF2Uper = par2

    def MPRLower(self,cMPR,cClosure,par1,par2,mode=2):
        if(cClosure < cMPR):
            print("WARNING cMPR shouldn't be bigger than cClosure")

        #here a bit of trickery cause eppler conventions are weird
        _N = self.ZeroN
        self.ZeroN = self.N / 2
        self.MPRLowNue = self.CtoNUpper(cMPR)
        self.cMPRlow = cMPR #for convinience
        self.CloseLowNue = self.CtoNUpper(cClosure)
        self.ZeroN = _N

        self.MPRModeLow = mode
        self.MPRF1Low = par1
        self.MPRF2Low = par2

    def writeMPR(self,mode,KR,Ktol=0) -> int:
        #we accept strings as input for mode but only 4 basis are avaliable
        #however all can be specified as integer
        if(isinstance(mode,str)): 
            if(mode == "upper"):
                mode = 1
            elif(mode == "lower"):
                mode = 2
            elif(mode == "mixed"):
                mode = 3
            elif(mode == "non"):
                mode = 0
            else:
                print("WARNING couldnt match mode to implemented type")
                return 1
        #write stuff to file
        File = self.File

        if(self.UseRamp != 0):  #the ramp hase been activated
            File.write("RAMP       ")
            #calculate delta nues
            dU = round(self.CtoNUpper(self.cMPRup - self.dc_Front_Up) - self.MPRUpNue,1)
            File.write(str(dU)+" ")
            dU = round(self.MPRUpNue - self.CtoNUpper(self.cMPRup + self.dc_Back_Up),1)
            File.write(str(dU)+" ")
            dU = round(self.CtoNLower(self.cMPRlow) - self.CtoNLower(self.cMPRlow - self.dc_Front_Low),1)
            File.write(str(dU)+" ")
            dU = round(self.CtoNLower(self.cMPRlow + self.dc_Back_Low) - self.CtoNLower(self.cMPRlow),1)
            File.write(str(dU)+"\n")
        
        File.write("TRA2" + "8" + " " + str(self.Number).zfill(4) + " ") #begin Tra2Line
        File.write(str(round(self.CloseUpNue,1)) + " " + str(round(self.MPRUpNue,1)) + " " + str(round(self.MPRModeUper,1)) + " " + str(round(self.MPRF1Uper,2)) + " " + str(round(self.MPRF2Uper,2)) + " " ) #Upper Side
        File.write(str(round(self.CloseLowNue,1)) + " " + str(round(self.MPRLowNue,1)) + " " + str(round(self.MPRModeLow,1)) + " " + str(round(self.MPRF1Low,2)) + " " + str(round(self.MPRF2Low,2)) + " " ) #lower side
        File.write(str(round(mode,0)) + " " + str(round(KR,1)) + " " + str(round(Ktol,1)) + "\n") #lower side
        
        return 0

    def RampUpper(self,dc_front,dc_back) -> int:
        self.UseRamp = 1
        self.dc_Front_Up = dc_front
        self.dc_Back_Up = dc_back
        return 0

    def RampLower(self,dc_front,dc_back):
        self.UseRamp = 1
        self.dc_Front_Low = dc_front
        self.dc_Back_Low = dc_back
        return 0

    def writeAlpha(self,alpha): #prints the ALFA command
        File = self.File
        File.write("ALFA     1 ")
        for i in range(0,len(alpha) - 1):
            File.write(str(round(alpha[i],1)) + " ")
        File.write(str(round(alpha[-1],1)) + "\n")

    def inviscidCalc_Custom(self,alphas) -> int:   #invicid calculation with array as alpha input
        self.writeAlpha(alphas)
        self.File.write("DIAG\n")
        return 0

    def invisidCalc(self,lowAlpha,upAlpha,n):   #invicid calc with range and number of points as input
        return self.inviscidCalc_Custom(onp.linspace(lowAlpha,upAlpha,num=n))

    def inviscidCalc_byIncrements(self,a0,da,n): #creates n points starting at a0 increasing by da
        x = onp.zeros(n)
        for i in range(0,n):
            x[i] = a0 + da * i
        return self.inviscidCalc_Custom(x)

    def naturalTransitionCalc(self,RE): #BL calculation with natural transition
        i = self.index_RE
        self.index_RE = i + 1
        self.REs[i] = RE #RE numbers for visous calculations
        self.index_transMode[i] = 3.01
        return 0

    def forcedTransitionCalc(self,RE,cTransUpper,cTransLower): #BL calculation with forced transition
        i = self.index_RE
        self.index_RE = i + 1
        self.REs[i] = RE #RE numbers for visous calculations
        self.index_transMode[i] = 1.01
        self.transTop[i] = cTransUpper
        self.transBottom[i] = cTransLower
        return 0

    def visousCalc(self,amin,amax,n):
        File = self.File
        self.writeAlpha(onp.linspace(amin,amax,num=n))
        File.write("RE  14     ")
        for i in range(0,5):
            if(self.index_transMode[i] == 0):
                File.write("0 " + str(round(self.REs[i],0)) + " ")
            else:
                File.write(str(self.index_transMode[i]) + " " + str(round(self.REs[i],0)) + " ")
        for i in range(0,4):
            File.write(str(round(self.transTop[i])) + " " + str(round(self.transBottom[i])) + " " )
        File.write(str(round(self.transTop[-1])) + " " + str(round(self.transBottom[-1])) + "\n" )
        File.write("CDCL\n")    
    
    def CloseFile(self) -> int: #Finalize the File
        self.File.write("ENDE")
        self.File.close()
        return 0