import sys
import numpy                      as np
import nkUtilities.LoadConfig     as lcf
import nkUtilities.plot1D         as pl1
import nkUtilities.configSettings as cfs


# ========================================================= #
# ===  display                                          === #
# ========================================================= #
def display( itpFile=None, refFile=None, pngFile=None, config=None ):
    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    if ( itpFile is None ): itpFile = "dat/xItp.dat"
    if ( refFile is None ): refFile = "dat/xRef.dat"
    if ( pngFile is None ): pngFile = "png/check__1d.png"
    if ( config  is None ): config  = lcf.LoadConfig()

    # ------------------------------------------------- #
    # --- [2] Fetch Data                            --- #
    # ------------------------------------------------- #
    with open( itpFile, "r" ) as f:
        iData = np.loadtxt( f )
    xItp =  iData[1:-1,0]
    yItp  = iData[1:-1,1] # edge 2 points is not good
    xAns  = iData[:   ,0]
    yAns  = iData[:   ,2]
    with open( refFile, "r" ) as f:
        rData = np.loadtxt( f )
    xRef  = rData[:,0]
    yRef  = rData[:,1]
    
    # ------------------------------------------------- #
    # --- [3] config Settings                       --- #
    # ------------------------------------------------- #
    cfs.configSettings( configType="plot1D_def", config=config )
    config["xTitle"]         = "X"
    config["yTitle"]         = "Y"
    config["plt_xAutoRange"] = True
    config["plt_yAutoRange"] = True
    config["plt_xRange"]     = [-5.0,+5.0]
    config["plt_yRange"]     = [-5.0,+5.0]
    config["plt_linewidth"]  = 1.0
    config["xMajor_Nticks"]  = 5
    config["yMajor_Nticks"]  = 5
    markersize               = 0.8

    # ------------------------------------------------- #
    # --- [4] plot Figure                           --- #
    # ------------------------------------------------- #
    fig = pl1.plot1D( config=config, pngFile=pngFile )
    fig.add__plot( xAxis=xItp, yAxis=yItp, label="Interpolation" )
    fig.add__plot( xAxis=xAns, yAxis=yAns, label="Analytical"    )
    fig.add__plot( xAxis=xRef, yAxis=yRef, label="Refference"   , \
                   linewidth=0.0, marker="o", markersize=markersize )
    fig.add__legend()
    fig.set__axis()
    fig.save__figure()


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    display()
