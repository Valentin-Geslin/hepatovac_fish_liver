def imageData = getCurrentImageData()
def name = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
double downsample = 1.0 
def labelServer = new LabeledImageServer.Builder(imageData)
    .backgroundLabel(0, ColorTools.WHITE)
    .downsample(downsample)   
    .addLabel('TISSUE', 1)  
    .multichannelOutput(false)
    .build()
annotations = getAnnotationObjects()
print annotations
removeObjects(annotations, true)
annotations.each{
    if(it.getPathClass() == null){return}
    addObject(it)
    className = it.getPathClass().toString()
    pathOutput = buildFilePath("/home/valentin/Desktop/hepatovac/tiles/", name, className)
    mkdirs(pathOutput)
    new TileExporter(imageData)
        .downsample(downsample)     
        .imageExtension('.png')     
        .tileSize(1600)      
        .labeledServer(labelServer)
        .annotatedTilesOnly(true)
        .overlap(0)          
        .writeTiles(pathOutput)
    removeObject(it, true)
    print 'Done!'
}
addObjects(annotations)
fireHierarchyUpdate()
