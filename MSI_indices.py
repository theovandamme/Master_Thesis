# Possible indices to use

def NDVI(image):
    # Calculate Normalized Difference Vegetation Index (NDVI)
    NDVI = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    return NDVI


def kNDVI(image):
    # Calculate Kernel NDVI (kNDVI)
    ndvi = NDVI(image)
    kNDVI = ndvi.pow(2).tanh().rename('kNDVI')
    return kNDVI

def NDWI(image):
    # Calculate Normalized Difference Water Index (NDWI)
    NDWI = image.normalizedDifference(['B3', 'B8']).rename('NDWI')
    return NDWI

def MNDWI(image):
    # Calculate Modified Normalized Difference Water Index (MNDWI)
    MNDWI = image.normalizedDifference(['B3', 'B11']).rename('MNDWI')
    return MNDWI

def RBR(image):
    # Calculate Red-Band Ratio (RBR)
    Red = image.select('B4')
    Green = image.select('B3')
    Blue = image.select('B2')
    RBR = (Red.divide(Red.add(Green).add(Blue))).rename('RBR')
    return RBR

def GBR(image):
    # Calculate Green-Band Ratio (GBR)
    Red = image.select('B4')
    Green = image.select('B3')
    Blue = image.select('B2')
    GBR = (Green.divide(Red.add(Green).add(Blue))).rename('GBR')
    return GBR

def EGI(image):
    # Calculate Excess Green Index (EGI)
    Red = image.select('B4')
    Green = image.select('B3')
    Blue = image.select('B2')
    EGI = (Green.multiply(ee.Number(2))).subtract(Red).subtract(Blue).rename('EGI')
    return EGI

def GRVI(image):
    # Calculate Green-Red Vegetation Index (GRVI)
    Red = image.select('B4')
    Green = image.select('B3')
    GRVI = (Green.subtract(Red)).divide(Green.add(Red)).rename('GRVI')
    return GRVI

def NDBRBI(image):
    # Calculate Sd Difference Bareness Index (NDBRBI)
    Red = image.select('B4')
    Blue = image.select('B2')
    NDBRBI = (Blue.subtract(Red)).divide(Blue.add(Red)).rename('NDBRBI')
    return NDBRBI

def SAVI(image):
    # Calculate Soil-Adjusted Vegetation Index (SAVI)
    SAVI = image.expression(
    '((NIR - RED) / (RED + NIR + L))*(1+L)',
    {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
        'L': ee.Number(0.5),
    },
    ).rename('SAVI')
    return SAVI

def TSAVI(image):
    # Calculate Transformed Soil-Adjusted Vegetation Index (TSAVI)
    TSAVI = image.expression(
    'a * ((NIR - a*RED - b) / (RED + a * NIR - a * b))',
    {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
        'a': ee.Number(0.33),
        'b': ee.Number(0.1),
    },
    ).rename('TSAVI')
    return TSAVI

def MSI(image):
    # Calculate Moisture Stress Index (MSI)
    NIR = image.select('B8')
    SWIR = image.select('B11')
    MSI = SWIR.divide(NIR).rename('MSI')
    return MSI

def LSWI(image):
    # Calculate Land Surface Water Index (LSWI)
    LSWI = image.normalizedDifference(['B8', 'B11']).rename('LSWI')
    return LSWI

def EVI(image):
    # Calculate Enhanced Vegetation Index (EVI)
    EVI = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RedEdge - 7.5 * BLUE + 1))',
    {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
        'BLUE': image.select('B2'),
        'RedEdge': image.select('B6')
    },
    ).rename('EVI')
    return EVI

def HSV(image):
    # Convert RGB image to HSV color space and extract the Hue component
    RGB = image.select(['B12','B8','B4'])
    HSV = RGB.rgbToHsv()
    H = HSV.select('hue').rename('HSV')
    return H

def HSV_1(image):
    # Convert image to HSV color space with custom calculations for H, S, and V components
    bands = ['B11', 'B8', 'B4']
    sentinel_image = image.select(bands)

    # Calculate V and S values
    V = sentinel_image.reduce(ee.Reducer.max())
    S = V.subtract(sentinel_image.reduce(ee.Reducer.min()))

    # Calculate H values using expressions
    H_R = sentinel_image.expression(
        '((nir - red) / S) * 60 + 60',  # Calculate H_R
        {'nir': sentinel_image.select('B8'),
         'red': sentinel_image.select('B4'),
         'S': S}
    )
    H_G = sentinel_image.expression(
        '((red - swir) / S) * 60 + 120',  # Calculate H_G
        {'red': sentinel_image.select('B4'),
         'swir': sentinel_image.select('B11'),
         'S': S}
    )
    H_B = sentinel_image.expression(
        '((swir - nir) / S) * 60 + 240',  # Calculate H_B
        {'swir': sentinel_image.select('B11'),
         'nir': sentinel_image.select('B8'),
         'S': S}
    )

    # Calculate H based on V and conditions
    H = V.where(V.eq(sentinel_image.select('B11')), H_R) \
        .where(V.eq(sentinel_image.select('B8')), H_G) \
        .where(V.eq(sentinel_image.select('B4')), H_B) \
        .where(V.eq(sentinel_image.reduce(ee.Reducer.min())), 0) \
        .rename('HSV_1')

    # Combine H, V, and S into an image
    HSV_1 = ee.Image([H, V, S]).rename(['HSV_1', 'V', 'S']).toUint16()
    return HSV_1

def MSAVI2(image):

    MSAVI2 = image.expression(
    '0.5*(2* NIR+1-sqrt((2*NIR+1)*(2*NIR+1)-8*(NIR-RED)))',
    {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
    },
    ).rename('MSAVI2')
    return MSAVI2

def TCT_Brightness(image):
  #https://doi.org/10.3390/rs10121862
  bands_to_resample = image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12'])
  # Resample all selected bands to 10m resolution
  resampled_image = bands_to_resample.resample('bilinear').reproject(
    crs=bands_to_resample.projection(),
    scale=10
  )
  result_image = image.addBands(resampled_image, overwrite=True)
  TCT_brightness = result_image.expression(
      '(0.0356*B1)+(0.0822*B2)+(0.1360*B3)+(0.2611*B4)+(0.2964*B5)+(0.338*B6)+(0.3877*B7)+(0.3895*B8)+(0.4750*B8A)+(0.949*B9)+(0.0009*B10)+(0.3882*B11)+(0.1366*B12)',
      {'B1' : image.select('B1'),
          'B2' : image.select('B2'),
          'B3' : image.select('B3'),
          'B4' : image.select('B4'),
          'B5' : image.select('B5'),
          'B6' : image.select('B6'),
          'B7' : image.select('B7'),
          'B8' : image.select('B8'),
          'B8A' : image.select('B8A'),
          'B9' : image.select('B9'),
          'B10' : image.select('B10'),
          'B11' : image.select('B11'),
          'B12' : image.select('B12'),},
  ).rename('TCT_Brightness')
  return TCT_brightness

def TCT_Greeness(image):
  #https://doi.org/10.3390/rs10121862
  bands_to_resample = image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12'])
  # Resample all selected bands to 10m resolution
  resampled_image = bands_to_resample.resample('bilinear').reproject(
    crs=bands_to_resample.projection(),
    scale=10
  )
  result_image = image.addBands(resampled_image, overwrite=True)
  TCT_Greeness = result_image.expression(
      '(*B1)+(B*B2)+(C*B3)+(D*B4)+(E*B5)+(F*B6)+(G*B7)+(H*B8)+(I*B8A)+(J*B9)+(K*B10)+(L*B11)+(M*B12)',
      {
          'B1' : image.select('B1'),
          'B2' : image.select('B2'),
          'B3' : image.select('B3'),
          'B4' : image.select('B4'),
          'B5' : image.select('B5'),
          'B6' : image.select('B6'),
          'B7' : image.select('B7'),
          'B8' : image.select('B8'),
          'B8A' : image.select('B8A'),
          'B9' : image.select('B9'),
          'B10' : image.select('B10'),
          'B11' : image.select('B11'),
          'B12' : image.select('B12'),

      },
  ).rename('TCT_Greeness')
  return TCT_Greeness
def TCT_Wetness(image):
  #https://doi.org/10.3390/rs10121862
  TCT_Wetness = image.expression(
      '(A*B1)+(B*B2)+(C*B3)+(D*B4)+(E*B5)+(F*B6)+(G*B7)+(H*B8)+(I*B8A)+(J*B9)+(K*B10)+(L*B11)+(M*B12)',
      {
          'B1' : image.select('B1'),
          'B2' : image.select('B2'),
          'B3' : image.select('B3'),
          'B4' : image.select('B4'),
          'B5' : image.select('B5'),
          'B6' : image.select('B6'),
          'B7' : image.select('B7'),
          'B8' : image.select('B8'),
          'B8A' : image.select('B8A'),
          'B9' : image.select('B9'),
          'B10' : image.select('B10'),
          'B11' : image.select('B11'),
          'B12' : image.select('B12'),
      },
  )

def TCT_Brightness_1(image):
  #https://doi.org/10.3390/rs10121862
  bands_to_resample = image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12'])
  # Resample all selected bands to 10m resolution
  resampled_image = bands_to_resample.resample('bilinear').reproject(
    crs=bands_to_resample.projection(),
    scale=10
  )
  result_image = image.addBands(resampled_image, overwrite=True)
  TCT_brightness_1 = result_image.expression(
      '(0.3037*B2)+(0.2793*B3)+(0.4743*B4)+(0.5585*B8)+(0.5082*B10)+(0.1863*B11)',
      {'B1' : image.select('B1'),
          'B2' : image.select('B2'),
          'B3' : image.select('B3'),
          'B4' : image.select('B4'),
          'B5' : image.select('B5'),
          'B6' : image.select('B6'),
          'B7' : image.select('B7'),
          'B8' : image.select('B8'),
          'B8A' : image.select('B8A'),
          'B9' : image.select('B9'),
          'B10' : image.select('B10'),
          'B11' : image.select('B11'),
          'B12' : image.select('B12'),},
  ).rename('TCT_Brightness_1')
  return TCT_brightness_1

def TCT_Greeness_1(image):
  #https://doi.org/10.3390/rs10121862
  bands_to_resample = image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12'])
  # Resample all selected bands to 10m resolution
  resampled_image = bands_to_resample.resample('bilinear').reproject(
    crs=bands_to_resample.projection(),
    scale=10
  )
  result_image = image.addBands(resampled_image, overwrite=True)
  TCT_brightness_1 = result_image.expression(
      '((-0.2848)*B2)+((-0.2435)*B3)+((-0.5436)*B4)+(0.7243*B8)+((-0.1800)*B11)+(0.0840*B12)',
      {'B1' : image.select('B1'),
          'B2' : image.select('B2'),
          'B3' : image.select('B3'),
          'B4' : image.select('B4'),
          'B5' : image.select('B5'),
          'B6' : image.select('B6'),
          'B7' : image.select('B7'),
          'B8' : image.select('B8'),
          'B8A' : image.select('B8A'),
          'B9' : image.select('B9'),
          'B10' : image.select('B10'),
          'B11' : image.select('B11'),
          'B12' : image.select('B12'),},
  ).rename('TCT_Brightness_1')
  return TCT_brightness_1
