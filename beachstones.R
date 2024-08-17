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

NAME="stone2"  # stone image filename
img=readPNG(paste0(NAME, ".png"))


# 1. CALCULATE 1px-WIDTH BORDER OF STONE
border=matrixfilter(img)
border[border==1]=0  # empty shape
border[border>0]=1  # set all detected borders to 1
border[img==0]=0  # drop borders not belonging to shape
writePNG(border, paste0(NAME, "_border.png"))


# 2. CALCULATE MAX DIAMETER
ind=which(border==1, arr.ind=TRUE)  # coords of pixels in border
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
x0=ind[p0,1]; y0=ind[p0,2]  # coords of stone max diameter
x1=ind[p1,1]; y1=ind[p1,2]
border=DrawLine(border, x0, y0, x1, y1, inc=FALSE, val=0.5)  # draw max diameter
writePNG(border, paste0(NAME, "_border_maxdiameter.png"))


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


# 4. CALCULATE TANGENTIAL POINTS
INI0=which(indorder[,1]==x0 & indorder[,2]==y0)
INI1=which(indorder[,1]==x1 & indorder[,2]==y1)
END0=INI1-1
END1=nrow(indorder)

m=(x0-x1)/(y1-y0)  # slope m=(alpha)
anglem=atan(m)*180/pi  # m in degrees

GAP=20
angulospre=c()
angulospos=c()
angulosmed=c()
mindist=9999999999
for (i in (INI0+GAP):(END0-GAP)) {
    xpre=indorder[i-GAP,1]
    ypre=indorder[i-GAP,2]
    
    xc=indorder[i,1]
    yc=indorder[i,2]
    
    xpos=indorder[i+GAP,1]
    ypos=indorder[i+GAP,2]
    
    tanalphapre=(xpre-xc)/(yc-ypre)
    tanalphapos=(xc-xpos)/(ypos-yc)
    
    anglepre=atan(tanalphapre)*180/pi
    anglepos=atan(tanalphapos)*180/pi    
    anglemed=(anglepre+anglepos)/2

    angulospre=c(angulospre, anglepre)    
    angulospos=c(angulospos, anglepos)
    angulosmed=c(angulosmed, anglemed)
    
    dist=abs(anglemed-anglem)
    if(dist < mindist) {
        mindist=dist
        pclosest=i
    }
    
}
plot(angulosmed, type='l')
lines(angulospre, col='red')
lines(angulospos, col='blue')

ma=-1/m
xa=indorder[pclosest,1]
ya=indorder[pclosest,2]

# Solve 2 equations linear system: A * k = b -> k = inv(A) * b
A=matrix(nrow=2, ncol=2)
A[1,]=c(m,  -1)
A[2,]=c(ma, -1)
b=as.matrix(c(m*x0-y0, ma*xa-ya))
xap=t(solve(A, b))  # equivalent to inv(A) * b = solve(A) %*% b
