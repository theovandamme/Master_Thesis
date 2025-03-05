import ee
import MSI_indices as indices

def extractIndex(img, index):
  index_img = indices.index(img)
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
