setImageType('BRIGHTFIELD_H_E');
setColorDeconvolutionStains('{"Name" : "H&E default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049", "Stain 2" : "Eosin", "Values 2" : "0.2159 0.8012 0.5581", "Background" : " 255 255 255"}');
runPlugin('qupath.imagej.detect.tissue.SimpleTissueDetection2', '{"threshold":220,"requestedPixelSizeMicrons":20,"minAreaMicrons":1.0E5,"maxHoleAreaMicrons":4500.0,"darkBackground":false,"smoothImage":false,"medianCleanup":false,"dilateBoundaries":false,"smoothCoordinates":true,"excludeOnBoundary":true,"singleAnnotation":true}')
resultingClass = getPathClass("TISSUE")
toChange = getAnnotationObjects().findAll{it.getPathClass() == null}
toChange.each{ it.setPathClass(resultingClass)}

import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathAnnotationObject

def project = getProject()
def imagesToExport = project.getImageList()
def separator = "\t"

def columnsToInclude = new String[]{""}
def exportType = PathAnnotationObject.class
def outputPath = "/home/valentin/Desktop/hepatovac/results/tissue_area_hepatovac.csv"
def outputFile = new File(outputPath)
def exporter  = new MeasurementExporter()
                  .imageList(imagesToExport)            // Images from which measurements will be exported
                  .separator(separator)                 // Character that separates values
                  .exportType(exportType)               // Type of objects to export
                  .exportMeasurements(outputFile)        // Start the export process
makeInverseAnnotation()
print "Done!"
