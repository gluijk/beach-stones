# Beach stones model
# www.overfitting.net
# https://www.overfitting.net/


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

for (names in 1:4) {
    NAME=paste0("stone", names)
    img=LoadBitmap(paste0(NAME, ".png"))
    
    # 1. CALCULATE 1px-WIDTH BORDER OF STONE
    border=matrixfilter(img)
    border[border==1]=0  # empty shape
    border[border>0]=1  # set all detected borders to 1
    border[img==0]=0  # drop borders not belonging to shape

    
    # 2. CALCULATE MAX DIAMETER
    ind=which(border==1, arr.ind=TRUE, useNames=FALSE)  # coords of pixels in border
    NUMPOINTS=nrow(ind)  # number of pixels in border
    maxdist=0
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
    border=DrawLine(border, x0, y0, x1, y1, inc=FALSE, val=0.5)  # draw max diameter
    border=DrawCircle(border, round(x0), round(y0),
                      RLARGE, inc=FALSE, val=0.5, thick=2)
    border=DrawCircle(border, round(x1), round(y1),
                      RLARGE, inc=FALSE, val=0.5, thick=2)
    
    
    # 3. ORDER BORDER PIXELS BY MIN NEIGHBOUR DISTANCE
    indorder=matrix(ind[p0,], nrow=1)  # ordered pixels starting from p0
    indcheck=ind[-p0,]  # remaining pixels to be checked
    
    # Recursively look for closest pixel to current checked pixel (i)
    for (i in 1:(nrow(indcheck)-1)) {
        NUMPOINTSCHECK=nrow(indcheck)  # number of border pixels still to check
        mindist=9999999999
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
    rm(indcheck)  # free tmp indcheck
    write.csv2(indorder, paste0(NAME, "_orderedpixels.csv"),
               row.names=FALSE)  # output ordered pixels
    
    
    # 4. CALCULATE 2 TANGENTIAL POINTS
    m=(y1-y0)/(x1-x0)  # slope m=tan(alpha)
    morth=-1/m         # slope orthogonal to m
    anglem=atan(m)*180/pi  # angle of m in degrees
    GAP=20  # number of pixels before and after tangential pixel checked
    
    NTAN=2  # number of tangential points to calculate
    START=unname(c(which(indorder[,1]==x0 & indorder[,2]==y0),
                   which(indorder[,1]==x1 & indorder[,2]==y1)))
    END=c(START[2]-1, nrow(indorder))
    
    for (n in 1:NTAN) {
        angulospre=c()
        angulospos=c()
        angulosmed=c()
        mindist=9999999999
        xrange=(START[n]+GAP):(END[n]-GAP)
        for (i in xrange) {
            # [[]] to unname all int values
            xpre=indorder[[i-GAP,1]]; ypre=indorder[[i-GAP,2]]
            xc=indorder[[i,1]];       yc=indorder[[i,2]]
            xpos=indorder[[i+GAP,1]]; ypos=indorder[[i+GAP,2]]
            
            tanalphapre=(yc-ypre)/(xc-xpre)
            tanalphapos=(ypos-yc)/(xpos-xc)
            
            anglepre=atan(tanalphapre)*180/pi
            anglepos=atan(tanalphapos)*180/pi    
            anglemed=(anglepre+anglepos)/2  # averate pre and pos calculated angles
        
            angulospre=c(angulospre, anglepre)    
            angulospos=c(angulospos, anglepos)
            angulosmed=c(angulosmed, anglemed)
            
            dist=abs(anglemed-anglem)
            if(dist < mindist) {
                mindist=dist
                pclosest=i
            }
            
        }
        png(paste0(NAME, "_tangentangles_",n,".png"), width=512, height=400)
        plot(xrange, angulosmed, type='l', ylim=c(-90,90), 
             main=paste0("Angles ", n ,"/2 of '", NAME,"'"),
             xlab="Contour pixel", ylab="Angle (º)", yaxt="n"
             )
        axis(2, at=seq(-90, 90, by=45), las=2)
        
        lines(xrange, angulospre, col='red')
        lines(xrange, angulospos, col='blue')
        abline(h=anglem, v=pclosest, lty=2)  # real angle of max diameter
        legend("topright", legend=c("Avg angle", "Pre angle", "Pos angle"),
               fill=c('black','red','blue'))
        dev.off()
    
        xorth=indorder[[pclosest,1]]
        yorth=indorder[[pclosest,2]]
        border=DrawCircle(border, round(xorth), round(yorth),
                          RSMALL, inc=FALSE, val=0.5)
        
        # Solve 2 equations linear system: A * x = b -> x = inv(A) * b
        # based on the equation of the line with a given slope (m and morth)
        # that passes through a point ((x0,y0) and (xorth,yorth)):
        # y-y0 = m*(x-x0)
        # y-yorth = morth*(x-xorth)
        A=matrix(nrow=2, ncol=2)
        A[1,]=c(m,     -1)
        A[2,]=c(morth, -1)
        b=as.matrix(c(m*x0-y0, morth*xorth-yorth))
        xorthp=t(solve(A, b))  # equivalent to inv(A) * b = solve(A) %*% b
    
        border=DrawLine(border, xorth, yorth, xorthp[1],
                        xorthp[2], inc=FALSE, val=0.5)  # draw max diameter
        border=DrawCircle(border, round(xorthp[1]), round(xorthp[2]),
                          RSMALL, inc=FALSE, val=0.5)
    }
    
    SaveBitmap(border, paste0(NAME, "_border_maxdiameter.png"))
}


