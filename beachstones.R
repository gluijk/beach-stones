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
    # Loop the kernel, vectorize the matrix
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

img=readPNG("stone2.png")


# Calculate 1px-width border of stone
border=matrixfilter(img)
border[border==1]=0  # empty shape
border[border>0]=1  # set all detected borders to 1
border[img==0]=0  # drop borders not belonging to shape

writePNG(border, "stoneborder2.png")
ind=which(border==1, arr.ind=TRUE)  # coords of pixels in border
NUMPOINTS=nrow(ind)  # number of pixels in border


# Order border pixels by min neighbour distance
indorder=matrix(ind[1,], nrow=1)  # ordered pixels
indcheck=ind[-1,]  # remaining pixels to be checked

for (i in 1:(nrow(indcheck)-1)) {
    NUMPOINTSCHECK=nrow(indcheck)  # number of border pixels to check
    mindist=9999999999
    indorderlast=indorder[nrow(indorder),]
    for (j in 1:NUMPOINTSCHECK) {
        dist=((indorderlast[1]-indcheck[j,1])^2 +
              (indorderlast[2]-indcheck[j,2])^2) ^ 0.5
        if(dist < mindist) {
            mindist=dist
            p1=j
        }
    }
    indorder=rbind(indorder, indcheck[p1,])  # add pixel to indorder
    indcheck=indcheck[-p1,]  # remove pixel from indcheck
}
indorder=rbind(indorder, indcheck)  # add last pixel to indorder
rm(indcheck)

write.csv2(indorder, "indorder2.csv")


# Calculate max diameter
maxdist=0
for (i in 1:(NUMPOINTS-1)) {
    for (j in (i+1):NUMPOINTS) {
        dist=((indorder[i,1]-indorder[j,1])^2 +
              (indorder[i,2]-indorder[j,2])^2) ^ 0.5
        if(dist > maxdist) {
            maxdist=dist
            p0=i
            p1=j
        }
    }
}
x0=indorder[p0,1]; y0=indorder[p0,2]  # coords of stone max diameter
x1=indorder[p1,1]; y1=indorder[p1,2]
border=DrawLine(border, x0, y0, x1, y1, inc=FALSE, val=0.5)  # draw max diameter
writePNG(border, "stoneborder_maxdiameter2.png")


# Calculate tangential points
tanalpha=(y1-y0)/(x1-x0)
angle=atan(tanalpha)*180/pi


