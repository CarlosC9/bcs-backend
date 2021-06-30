import unittest
import requests
import json

class RegionsAPITest(unittest.TestCase):
    url = "http://localhost:5000/api/geo/regions/"
    post_req =  json.dumps({"attributes": {"tags": ["test", "testing"]},
         "name": "newnew",
         "geometry": {
             "type": "FeatureCollection",
             "features": [
                 {
                     "type": "Feature",
                     "properties": {},
                     "geometry": {
                         "type": "Polygon",
                         "coordinates": [
                             [
                                 [
                                     -15.932922363281252,
                                     28.282614623837407
                                 ],
                                 [
                                     -15.92742919921875,
                                     28.28140526783416
                                 ],
                                 [
                                     -16.489105224609375,
                                     28.33943885710451
                                 ],
                                 [
                                     -16.446533203125,
                                     27.925260662705618
                                 ],
                                 [
                                     -15.932922363281252,
                                     28.282614623837407
                                 ]
                             ]
                         ]
                     }
                 }
             ]
         }
    })


    put_req_1 =  json.dumps({  "tags" : ["test","testing"],
                    "name" : "changed"})

    put_req_2 = json.dumps({"geometry": {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [
                                -15.932922363281252,
                                28.282614623837407
                            ],
                            [
                                -15.92742919921875,
                                28.28140526783416
                            ],
                            [
                                -16.489105224609375,
                                28.33943885710451
                            ],
                            [
                                -16.5,
                                27.925260662705618
                            ],
                            [
                                -15.932922363281252,
                                28.282614623837407
                            ]
                        ]
                    ]
                }
            }
        ]
    }
    })

    session = {'Cookie': 'session=242b93a0-ae2d-4961-86f3-2de20c9bb3a6'}

    def test_log_in(self):
        url = "http://localhost:5000/api/authn?user=test_user"

        response = requests.request("PUT", url)

        print(response.text)
        return response.cookies

    def test_create_region(self):
        url = "http://localhost:5000/api/geo/regions/"

        payload = self.post_req
        headers = self.session

        response = requests.request("POST", self.url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def test_get_region_by_id(self):
        import requests

        id = 1

        payload = {}
        headers = self.session

        response = requests.request("GET", f"{self.url}{id}", headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def test_get_regions(self):

        payload = {}
        headers = {
            'Cookie': 'session=242b93a0-ae2d-4961-86f3-2de20c9bb3a6'
        }

        response = requests.request("GET", self.url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def test_modify_region_attributes(self):

        id = 1

        payload = self.put_req_1
        headers = self.session

        response = requests.request("PUT", f"{self.url}{id}", headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None


    def test_modify_region_geometry(self):

        url = "http://localhost:5000/api/geo/regions/1"

        payload = self.put_req_2
        headers = {
            'Cookie': 'session=242b93a0-ae2d-4961-86f3-2de20c9bb3a6'
        }

        response = requests.request("PUT", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code == 200)

        return None

    def test_detele_region(self):

        url = "http://localhost:5000/api/geo/regions/1"

        payload = {}
        headers = {
            'Cookie': 'session=242b93a0-ae2d-4961-86f3-2de20c9bb3a6'
        }

        response = requests.request("DELETE", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code == 200)

        return None


class LayersAPITest(unittest.TestCase):
    session = {'Cookie': 'session=242b93a0-ae2d-4961-86f3-2de20c9bb3a6'}

    post_req = json.dumps({"data": {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "Temp": 20,
                    "ind": 5
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [
                                -15.450897216796873,
                                28.121952659013388
                            ],
                            [
                                -15.435104370117188,
                                28.078341712094886
                            ],
                            [
                                -15.412445068359377,
                                28.11014311115806
                            ],
                            [
                                -15.443344116210936,
                                28.16251913936131
                            ],
                            [
                                -15.450897216796873,
                                28.121952659013388
                            ]
                        ]
                    ]
                }
            },
            {
                "type": "Feature",
                "properties": {
                    "Temp": 30,
                    "ind": 6
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [
                                -15.618438720703123,
                                28.16584854109881
                            ],
                            [
                                -15.739288330078123,
                                28.10166365943778
                            ],
                            [
                                -15.65826416015625,
                                28.075615439518966
                            ],
                            [
                                -15.618438720703123,
                                28.16584854109881
                            ]
                        ]
                    ]
                }
            },
            {
                "type": "Feature",
                "properties": {
                    "Temp": 40,
                    "ind": 7
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [
                                -15.812759399414062,
                                27.880963078302393
                            ],
                            [
                                -15.605392456054686,
                                27.73945373313175
                            ],
                            [
                                -15.391845703125,
                                27.881570017022806
                            ],
                            [
                                -15.812759399414062,
                                27.880963078302393
                            ]
                        ]
                    ]
                }
            }
        ]
    },
        "name": "custom_data",
        "wks": "ngd"
    })

    def test_create_layer_by_vector_file(self):

        url = "http://localhost:5000/api/geo/layers/"
        payload = {'data': '{"name":"test"}'}
        files = {
            'file1': (
                'Plantas.shp', open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/plantae_canarias/Plantas.shp', 'rb'),
                'application/octet-stream'),
            'file2': (
                'Plantas.shx', open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/plantae_canarias/Plantas.shx', 'rb'),
                'application/octet-stream'),
            'file3': (
                'Plantas.dbf', open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/plantae_canarias/Plantas.dbf', 'rb'),
                'application/octet-stream'),
            'file4': (
                'Plantas.prj', open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/plantae_canarias/Plantas.prj', 'rb'),
                'application/octet-stream')
        }
        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload, files=files)

        print(response.text)

        self.assertEqual(response.status_code, 200)
        return  None

    def test_create_layer_by_raster_file(self):
        url = "http://localhost:5000/api/geo/layers/"

        payload = {'data': '{"name": "raster_test"}'}
        files = [
            ('file', ('RT_SIS.tif',
                      open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/gobcan_riesgomap_cartografia/RT_SIS.tif', 'rb'),
                      'image/tiff')),
            ('', ('RT_SIS.tfw',
                  open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/gobcan_riesgomap_cartografia/RT_SIS.tfw', 'rb'),
                  'application/octet-stream'))
        ]
        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload, files=files)

        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_create_layer_by_geojson_data(self):


        payload = self.post_req

        url = "http://localhost:5000/api/geo/layers/"

        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_create_convert_layer(self):
        """
        use a biota shape file as input (Plantas_canarias)
        and convert it into JSON for PDA study
        """

        url = "http://localhost:5000/api/geo/layers/"

        payload = {
            'data': '{"name":"plantas_test1","attributes": {"tags":["plantas","Canarias","Biota"]},"wks":"ngd","convert_to":"geojson"}'}
        files = [
            ('file1', (
            'Plantas.prj', open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/plantae_canarias/Plantas.prj', 'rb'),
            'application/octet-stream')),
            ('file2', (
            'Plantas.shx', open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/plantae_canarias/Plantas.shx', 'rb'),
            'application/octet-stream')),
            ('file3', (
            'Plantas.dbf', open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/plantae_canarias/Plantas.dbf', 'rb'),
            'application/octet-stream')),
            ('file4', (
            'Plantas.shp', open('/home/paula/Documentos/NEXTGENDEM/tests_GIS/plantae_canarias/Plantas.shp', 'rb'),
            'application/octet-stream'))
        ]
        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload, files=files)

        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_get_layer_by_id(self):

        id = 1

        url = f"http://localhost:5000/api/geo/layers/{id}"

        payload = {}
        headers = self.session

        response = requests.request("GET", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)
        return None

    def test_get_layers(self):

        url = "http://localhost:5000/api/geo/layers/"

        payload = {}
        headers = self.session

        response = requests.request("GET", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def test_create_temporal_view(self):
        """
        create a new view in geoserver from a existing layer (one or more) stored in postgis
        using a query.
        This layer will be temporal ans will be deleted any time that a new temporal layer is created
        @return:
        """

        sql = "SELECT geometry as geom, \"ind\" FROM layer_31 WHERE \"ind\" =7"
        key_col = "ind"
        url = f"http://localhost:5000/api/geo/layers/?filter={sql}&key_col={key_col}"

        payload = {}
        headers = {
            'Cookie': 'session=242b93a0-ae2d-4961-86f3-2de20c9bb3a6'
        }

        response = requests.request("GET", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def test_save_temportal_view(self):

        """
        it convert temporal views into vector layers in postgis and publish it in geoserver
        @return:
        """


        url = "http://localhost:5000/api/geo/layers/"

        payload = "{ \"layer_name\": \"tmpview\",\"name\": \"test_tmpview\", \"wks\":\"ngd\"}"
        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def test_create_user_view(self):
        layer_name = "user_view"
        sql = "SELECT geometry as geom, \"ind\" FROM layer_31 WHERE \"ind\" =7"
        key_col = "ind"
        url = f"http://localhost:5000/api/geo/layers/?filter={sql}&key_col={key_col}"

        payload = {}
        headers = {
            'Cookie': 'session=242b93a0-ae2d-4961-86f3-2de20c9bb3a6'
        }

        response = requests.request("GET", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def test_save_user_view(self):
        """
        it convert user views into vector layers in postgis and publish it in geoserver
        """
        url = "http://localhost:5000/api/geo/layers/"

        layer_name = "user_view"

        payload = "{ \"layer_name\": \"tmpview_3\",\"name\": \"test_tmpview\", \"wks\":\"ngd\"}"
        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return  None

    def test_delete_layer(self):
        """
        @param id:
        @return:
        """

        id = 1

        url = f"http://localhost:5000/api/geo/layers/{id}"

        payload = {}
        headers = self.session

        response = requests.request("DELETE", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

class StylesAPITest(unittest.TestCase):

    session = {'Cookie': 'session=242b93a0-ae2d-4961-86f3-2de20c9bb3a6'}


    def raster_style_ramp_by_value(self):
        """
        change color ramp in a already created style for a raster layer
        @return:
        """


        url = "http://localhost:5000/api/geo/styles/"

        payload  = json.dumps({"c_ramp_1" : {
            'label 1 value': '#ffff55',
            'label 2 value': '#505050',
            'label 3 value': '#404040',
            'label 4 value': '#333333'
            },
            "file" : "path/to/raster/file"
        })
        headers = self.session


        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None


        pass

    def raster_style_by_defined(self):
        """
        change color ramp in a already created style for a raster layer
        @return:
        """

        url = "http://localhost:5000/api/geo/styles/"

        # (https://matplotlib.org/3.3.0/tutorials/colors/colormaps.html)

        payload = json.dumps(
            {"c_ramp_1": "RdGy" ,
            "file": "path/to/raster/file"
        })

        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def raster_style_ramp_by_list(self):
        url = "http://localhost:5000/api/geo/styles/"

        payload = json.dumps(
            {"c_ramp_1": ["#ffffff", "#453422",  "#f0f0f0", "#aaaaaa"],
             "file": "path/to/raster/file"
             })

        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def vector_style_for_boundaries(self):
        """add a style for the boundary in this case, multypoligon"""

        url = "http://localhost:5000/api/geo/styles/"

        payload = json.dumps(
            {"color": "#ffffff",
             "style_mame": "new_style"
             })

        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def vector_style_for_features_data(self):
        """
        add a style for the boundary in this case, multypoligon
        """


        url = "http://localhost:5000/api/geo/styles/"

        payload = json.dumps(
            {"column": "isla", #mandatory
             "categorized_data": ["Gran Canaria","Tenerife"],
             "color_ramp" : ["#ffff55", "#505050'"]
             })

        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None

    def vector_style_for_features_number(self):
        """
        add a style for the boundary in this case, multypoligon
        """


        url = "http://localhost:5000/api/geo/styles/"

        payload = json.dumps(
            {"column": "RIQUEZA", #mandatory
             "color_ramp" : "RdGy"
             }
        )

        headers = self.session

        response = requests.request("POST", url, headers=headers, data=payload)

        print(response.text)

        self.assertEqual(response.status_code, 200)

        return None



if __name__ == '__main__':
    RegionsAPITest.test_create_region()
    RegionsAPITest.test_get_regions()
    RegionsAPITest.test_modify_region_geometry()
    RegionsAPITest.test_modify_region_attributes()
    RegionsAPITest.test_detele_region()


    LayersAPITest.test_create_layer_by_geojson_data()
    LayersAPITest.test_create_layer_by_vector_file()
    LayersAPITest.test_create_convert_layer()
    LayersAPITest.test_create_layer_by_raster_file()
    LayersAPITest.test_create_temporal_view()


