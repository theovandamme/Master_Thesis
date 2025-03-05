import ee
import MSI_indices as indices
import numpy as np

def extractIndex(img, index):
  print(f"Index passed: {'NDVI'}")
  index_function = getattr(indices, 'NDVI')
  index_img = index_function(img)
  band_name = index_img.bandNames().get(0)
  calc = index_img.select(ee.String(band_name)).reduceRegion(
        reducer=ee.Reducer.toList(),
        geometry=region,
        scale=10
    ).get(ee.String(band_name))
  return ee.Feature(None, {
        'index': calc,
        'date': img.date().format('YYYY-MM-dd')
    })

def convert_to_numpy(collection, index):
  print(index)
  #index_collection = collection.map(index)
  # Extract NDVI values and dates
  index_values = collection.map(lambda img: extractIndex(img, index))
  # Get the data as a list
  index_list = index_values.aggregate_array('index').getInfo()
  dates_list = index_values.aggregate_array('date').getInfo()

  # Convert to NumPy array
  index_array = np.array(index_list)
  dates_array = np.array(dates_list)

  index_array = index_array.T
  # Print results
  return index_array, dates_array
