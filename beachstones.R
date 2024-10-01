# Beach stones model
# www.overfitting.net
# https://www.overfitting.net/2024/08/modelo-geometrico-de-piedra-de-playa.html


library(png)

# Por Carlos Gil Bellosta
indices.drawline = function(x0, y0, x1, y1) {
    x0=round(x0)
    x1=round(x1)
    y0=round(y0)
    y1=round(y1)
    
    if (y0 == y1) return(cbind(x0:x1, y0)) # Recta de m=0 o un punto
    if (abs(x1 - x0) >= abs(y1 - y0)) { # Recta de 0 < |m| <= 1
        m = (y1 - y0) / (x1 - x0)
        cbind(x0:x1, round(y0 + m * ((x0:x1) - x0)))
    } else indices.drawline(y0, x0, y1, x1)[, 2:1]  # Recta de |m| > 1
    # Llamada traspuesta recursiva y traspuesta
}

DrawLine = function(img, x0, y0, x1, y1, inc=TRUE, val=1) {
    # Dibuja recta desde (x0,y0)-(x1,y1)
    # Por defecto método no destructivo y con valor=1
    indices=indices.drawline(x0, y0, x1, y1)
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}

DrawPoint = function(img, x0, y0, inc=TRUE, val=1) {
    # Dibuja punto en (x0,y0)
    # Por defecto método no destructivo y con valor=1
    img=DrawLine(img, x0, y0, x0, y0, inc, val)
    
    return(img)
}

DrawEllip = function(img, x0, y0, a, b, inc=TRUE, val=1, fill=FALSE, thick=1) {
    # Dibuja elipse de centro (x0,y0) y radios a y b
    # Por defecto método no destructivo, con valor=1 y sin relleno
    # Puede elegirse el grosor si no se rellena
    # Aquí no redondeamos para tener más precisión en la división
    if (fill) {
        indices=which( ((row(img)-x0)/a)^2 + ((col(img)-y0)/b)^2 < 1 )
    } else {
        indices=which( ((row(img)-x0)/(a+thick/2))^2 + ((col(img)-y0)/(b+thick/2))^2 <  1 &
                       ((row(img)-x0)/(a-thick/2))^2 + ((col(img)-y0)/(b-thick/2))^2 >= 1 )
    }
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}

DrawCircle = function(img, x0, y0, r, inc=TRUE, val=1, fill=FALSE, thick=1) {
    # Dibuja círculo de centro (x0,y0) y radio r
    # Por defecto método no destructivo, con valor=1 y sin relleno
    # Puede elegirse el grosor si no se rellena
    img=DrawEllip(img, x0, y0, r, r, inc, val, fill, thick)
    
    return(img)
}

LoadBitmap = function(name, chan=2) {
    # Lee bitmap en formato PNG
    # Si no es monocromo se carga el canal chan (por defecto G)
    require(png)
    img=readPNG(name)
    if (length(dim(img))>2) img=img[,,chan]
    
    return(t(img[nrow(img):1,]))
}

SaveBitmap = function(img, name, trunc=TRUE, gamma=1) {
    # Guarda bitmap en formato PNG
    # Solo si trunc=FALSE y la imagen excede de 1 se reescala a 1
    require(png)
    img[img<0]=0
    if (trunc) img[img>1]=1
    if (tolower(substr(name, nchar(name)-3, nchar(name))) != ".png") name=paste0(name,".png")
    writePNG(t(img[,ncol(img):1] / max(max(img),1))^(1/gamma), name)
}


# Matrix vectorized convolution by kernel
matrixfilter=function(img, kernel=matrix(c(0,0.25,0, 0.25,0,0.25, 0,0.25,0),
                                         nrow=3, ncol=3)) {
    KY=nrow(kernel); DIMY=nrow(img)
    KX=ncol(kernel); DIMX=ncol(img)
    
    imgoutpre=matrix(0, nrow=DIMY-KY+1, ncol=DIMX-KX+1)  # inner subset of img
    # Convolution: loop the kernel/vectorize the matrix
    for (j in 1:KY) {
        for (i in 1:KX) {
            if (kernel[j,i]) imgoutpre=imgoutpre +
                    img[(1+j-1):(j+DIMY-KY), (1+i-1):(i+DIMX-KX)]*kernel[j,i]
        }
    }
    
    imgout=img  # keep unfiltered img borders
    imgout[(1+(KY-1)/2):(DIMY-(KY-1)/2),  # overwrite filtered area
           (1+(KX-1)/2):(DIMX-(KX-1)/2)]=imgoutpre

    return (imgout)
}


####################################

RLARGE=7  # circles to plot radius
RSMALL=5
NSTONES=2

for (stone in 1:NSTONES) {  # process stones
    NAME=paste0("stone", stone)
    img=LoadBitmap(paste0(NAME, ".png"))


    # 1. CALCULATE 1px-WIDTH BORDER OF STONE
    border=matrixfilter(img)
    border[border==1]=0  # empty shape
    border[border>0]=1  # set all detected borders to 1
    border[img==0]=0  # drop borders not belonging to shape
    DIMX=nrow(border)
    DIMY=ncol(border)
    
    
    # 2. CALCULATE MAX DIAMETER
    ind=which(border==1, arr.ind=TRUE, useNames=FALSE)  # coords of pixels in border
    NUMPOINTS=nrow(ind)  # number of pixels in border
    maxdist=0  #arbitrarily low number
    for (i in 1:(NUMPOINTS-1)) {
        for (j in (i+1):NUMPOINTS) {
            dist=((ind[i,1]-ind[j,1])^2 +
                  (ind[i,2]-ind[j,2])^2) ^ 0.5
            if(dist > maxdist) {
                maxdist=dist
                p0=i
                p1=j
            }
        }
    }
    x0=ind[p0,1]  # coords of stone max diameter
    y0=ind[p0,2]
    x1=ind[p1,1]
    y1=ind[p1,2]
    onlylines=border*0
    onlylines=DrawLine(onlylines, x0, y0, x1, y1,
                       inc=FALSE, val=0.5)  # draw max diameter
    onlylines=DrawCircle(onlylines, round(x0), round(y0),
                         RLARGE, inc=FALSE, val=0.5, thick=2)
    onlylines=DrawCircle(onlylines, round(x1), round(y1),
                         RLARGE, inc=FALSE, val=0.5, thick=2)
    
    
    # 3. ORDER BORDER PIXELS BY MIN NEIGHBOUR DISTANCE
    indorder=matrix(ind[p0,], nrow=1)  # ordered pixels starting from p0
    indcheck=ind[-p0,]  # remaining pixels to be checked
    
    # Recursively look for closest pixel to current checked pixel (i)
    for (i in 1:(nrow(indcheck)-1)) {
        NUMPOINTSCHECK=nrow(indcheck)  # number of border pixels still to check
        mindist=9999999999  # arbitrarily high number
        indorderlast=indorder[nrow(indorder),]  # last pixel added to indorder
        for (j in 1:NUMPOINTSCHECK) {
            dist=((indorderlast[1]-indcheck[j,1])^2 +
                  (indorderlast[2]-indcheck[j,2])^2) ^ 0.5
            if(dist < mindist) {
                mindist=dist
                pclosest=j
            }
        }
        indorder=rbind(indorder, indcheck[pclosest,])  # add pixel to indorder
        indcheck=indcheck[-pclosest,]  # remove pixel from indcheck
    }
    indorder=rbind(indorder, indcheck)  # add last pixel to indorder
    write.csv2(indorder, paste0(NAME, "_orderedpixels.csv"),
               row.names=FALSE)  # output ordered pixels
    
    
    # 4. CALCULATE 2 TANGENTIAL POINTS AND DRAW 4 ELLIPSES QUARTERS
    m=(y1-y0)/(x1-x0)  # slope of max diameter (m=tan(alpha))
    morth=-1/m         # slope orthogonal to m
    anglemrad=atan(m)        # angle of m in radians
    anglem=anglemrad*180/pi  # angle of m in degrees
    GAP=40  # number of pixels before and after tangential pixel checked
    
    # Starting and ending pixel index to check on each side of the diameter
    START=unname(c(which(indorder[,1]==x0 & indorder[,2]==y0),
                   which(indorder[,1]==x1 & indorder[,2]==y1)))
    END=c(START[2]-1, nrow(indorder))
    
    model=onlylines*0
    for (n in 1:2) {  # 2 tangential points to calculate
        angulospre=c()
        angulospos=c()
        angulosmed=c()
        mindist=9999999999  #arbitrarily high number
        xrange=(START[n]+GAP):(END[n]-GAP)
        for (i in xrange) {
            # all [[]] are to unname all int values
            xpre=indorder[[i-GAP,1]]; ypre=indorder[[i-GAP,2]]
            xc=indorder[[i,1]];       yc=indorder[[i,2]]
            xpos=indorder[[i+GAP,1]]; ypos=indorder[[i+GAP,2]]
            
            tanalphapre=(yc-ypre)/(xc-xpre)
            tanalphapos=(ypos-yc)/(xpos-xc)
    
            # Averaging pre and pos angles is much more
            # convenient than averaging pre and pos slopes            
            anglepre=atan(tanalphapre)*180/pi
            anglepos=atan(tanalphapos)*180/pi    
            anglemed=(anglepre+anglepos)/2
    
            angulospre=c(angulospre, anglepre)    
            angulospos=c(angulospos, anglepos)
            angulosmed=c(angulosmed, anglemed)
            
            dist=abs(anglemed-anglem)
            if(dist < mindist) {  # pixel with tangent angle closest to anglem
                mindist=dist
                pclosest=i
            }
            
        }
    
        xorth=indorder[[pclosest,1]]
        yorth=indorder[[pclosest,2]]
        onlylines=DrawCircle(onlylines, round(xorth), round(yorth),
                             RSMALL, inc=FALSE, val=0.5)
        
        # Draw angle for every contour pixel
        png(paste0(NAME, "_tangentangles_",n,".png"), width=512, height=400)
        plot(xrange, angulosmed, type='l', ylim=c(-90,90), 
             main=paste0("Angles ", n ,"/2 of '", NAME,"'"),
             xlab="Contour pixel", ylab="Angle (º)", yaxt="n"
             )
        axis(2, at=seq(-90, 90, by=45), las=2)
        lines(xrange, angulospre, col='red')
        lines(xrange, angulospos, col='blue')
        abline(h=anglem, v=pclosest, lty=2, col='gray')  # real angle of max diameter
        legend("topright", legend=c("Avg angle", "Pre angle", "Pos angle"),
               fill=c('black','red','blue'))
        dev.off()
    
        # Solve 2 equations linear system: A * k = b -> k = inv(A) * b
        # based on the equation of the line with a given slope
        # that passes through a point:
        #   y-y0    = m     * (x-x0)
        #   y-yorth = morth * (x-xorth)
        A=matrix(nrow=2, ncol=2)
        A[1,]=c(m,     -1)
        A[2,]=c(morth, -1)
        b=as.matrix(c(m*x0-y0, morth*xorth-yorth))
        k=t(solve(A, b))  # equivalent to inv(A) * b = solve(A) %*% b
        xorthc=k[1]
        yorthc=k[2]
    
        onlylines=DrawLine(onlylines, xorth, yorth, xorthc, yorthc,
                           inc=FALSE, val=0.5)  # draw max diameter
        onlylines=DrawCircle(onlylines, round(xorthc), round(yorthc),
                             RSMALL, inc=FALSE, val=0.5)
        
        # Draw 2 ellipses defined by (xorthc, yorthc) and (x0,y0) and (x1,y1)
        RESOL=100  # number of points to define each of the 1/4 ellipses
        # Ellipse B axis
        axisB=((xorthc-xorth)^2+(yorthc-yorth)^2)^0.5
        ellipX=c(x0,x1)
        ellipY=c(y0,y1)
        for (ellip in 1:2) {  # 2 partial ellipses to plot for each tan point
            # Ellipse A axis
            axisA=((xorthc-ellipX[ellip])^2+(yorthc-ellipY[ellip])^2)^0.5
            anglerange=seq(0, pi/2, pi/2/RESOL)
            
            # Ad-hoc taylored 4 sections of ellipses
            # (room for improvement)
            if (n==1 & ellip==2) anglerange=anglerange+pi/2
            if (n==2 & ellip==1) anglerange=anglerange-pi/2
            if (n==2 & ellip==2) anglerange=anglerange+pi
    
            xpre=axisA*cos(anglerange[1])
            ypre=axisB*sin(anglerange[1])
            xplotpre=xorthc+xpre*cos(-anglemrad)+ypre*sin(-anglemrad)
            yplotpre=yorthc-xpre*sin(-anglemrad)+ypre*cos(-anglemrad)
            for (angle in anglerange) {
                xpre=axisA*cos(angle)
                ypre=axisB*sin(angle)
                # Ellipse 2D rotation by m (not sure why anglemrad needs -)
                xplot=xorthc+xpre*cos(-anglemrad)+ypre*sin(-anglemrad)
                yplot=yorthc-xpre*sin(-anglemrad)+ypre*cos(-anglemrad)
                if (xplot>=0 & xplot<=DIMX & yplot>=0 & yplot<=DIMY)
                    model=DrawLine(model, xplotpre,yplotpre, xplot,yplot,
                                   inc=FALSE, val=1)
                xplotpre=xplot
                yplotpre=yplot
            }
        }
    }
    
    SaveBitmap(model, paste0(NAME, "_model_only.png"))
    SaveBitmap(border, paste0(NAME, "_border_only.png"))
    SaveBitmap(onlylines/max(onlylines), paste0(NAME, "_lines_only.png"))
    
    border[onlylines > 0]=onlylines[onlylines > 0]
    model[onlylines > 0]=onlylines[onlylines > 0]
    SaveBitmap(border, paste0(NAME, "_border_lines.png"))    
    SaveBitmap(model, paste0(NAME, "_model_lines.png"))
}


# 5. CHECK MODEL FITTING

# First manually create solid versions of each model (fill) -> _model_solid.png
for (stone in 1:NSTONES) {  # process stones
    NAME=paste0("stone", stone)

    # Solid versions of stone and stone model
    bordersolid=LoadBitmap(paste0(NAME, ".png"))
    modelsolid=LoadBitmap(paste0(NAME, "_model_solid.png"))

    intersec=length(which(bordersolid==1 & modelsolid==1))  # common pixels
    # A: % of the model shape that falls into the real stone
    A=intersec/sum(modelsolid)
    # B: % of the real stone that falls into the model shape
    B=intersec/sum(bordersolid)
    Accuracy=(A * B)^0.5
    
    print(paste0("Performance of '", NAME, "' model: ",
                 "A=", round(A*100,2), "%, ",
                 "B=", round(B*100,2), "%, ",
                 "Accuracy=(A*B)^0.5=", round(Accuracy*100,2), "%"))
}
