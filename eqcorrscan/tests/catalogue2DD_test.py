"""
Functions to test the functions within the eqcorrscan.utils.catalogue2DD.py \
submodule.  Uses test data distributed with the EQcorrscan pacakge.
"""
import unittest


class TestCatalogueMethods(unittest.TestCase):
    def test_rounding(self):
        """Simple test to test that _cc_round gives correct result.
        """
        from eqcorrscan.utils.catalogue2DD import _cc_round
        import numpy as np
        # Use an irrational number and round it to various decimal places
        test_no = np.pi
        self.assertEqual(_cc_round(test_no, 0), '3')
        self.assertEqual(_cc_round(test_no, 1), '3.1')
        self.assertEqual(_cc_round(test_no, 2), '3.14')
        self.assertEqual(_cc_round(test_no, 3), '3.142')
        self.assertEqual(_cc_round(test_no, 4), '3.1416')
        self.assertEqual(_cc_round(test_no, 5), '3.14159')
        self.assertEqual(_cc_round(test_no, 6), '3.141593')
        self.assertEqual(_cc_round(test_no, 7), '3.1415927')
        self.assertEqual(_cc_round(test_no, 8), '3.14159265')
        self.assertEqual(_cc_round(test_no, 9), '3.141592654')
        self.assertEqual(_cc_round(test_no, 10), '3.1415926536')
        self.assertEqual(_cc_round(test_no, 11), '3.14159265359')
        self.assertEqual(_cc_round(test_no, 12), '3.141592653590')
        self.assertEqual(_cc_round(test_no, 13), '3.1415926535898')
        self.assertEqual(_cc_round(test_no, 14), '3.14159265358979')
        self.assertEqual(_cc_round(test_no, 15), '3.141592653589793')

    def test_weight_averaging(self):
        """Simple function to test _av_weight returns the correct weights.
        """
        from eqcorrscan.utils.catalogue2DD import _av_weight
        self.assertEqual(_av_weight('0', '0'), '1.0000')
        self.assertEqual(_av_weight('0', '1'), '0.8750')
        self.assertEqual(_av_weight('0', '2'), '0.7500')
        self.assertEqual(_av_weight('0', '3'), '0.6250')
        self.assertEqual(_av_weight('0', '4'), '0.5000')
        self.assertEqual(_av_weight('0', '9'), '0.5000')
        self.assertEqual(_av_weight('1', '0'), '0.8750')
        self.assertEqual(_av_weight('1', '1'), '0.7500')
        self.assertEqual(_av_weight('1', '2'), '0.6250')
        self.assertEqual(_av_weight('1', '3'), '0.5000')
        self.assertEqual(_av_weight('1', '4'), '0.3750')
        self.assertEqual(_av_weight('1', '9'), '0.3750')
        self.assertEqual(_av_weight('2', '0'), '0.7500')
        self.assertEqual(_av_weight('2', '1'), '0.6250')
        self.assertEqual(_av_weight('2', '2'), '0.5000')
        self.assertEqual(_av_weight('2', '3'), '0.3750')
        self.assertEqual(_av_weight('2', '4'), '0.2500')
        self.assertEqual(_av_weight('2', '9'), '0.2500')
        self.assertEqual(_av_weight('3', '0'), '0.6250')
        self.assertEqual(_av_weight('3', '1'), '0.5000')
        self.assertEqual(_av_weight('3', '2'), '0.3750')
        self.assertEqual(_av_weight('3', '3'), '0.2500')
        self.assertEqual(_av_weight('3', '4'), '0.1250')
        self.assertEqual(_av_weight('3', '9'), '0.1250')
        self.assertEqual(_av_weight('4', '0'), '0.5000')
        self.assertEqual(_av_weight('4', '1'), '0.3750')
        self.assertEqual(_av_weight('4', '2'), '0.2500')
        self.assertEqual(_av_weight('4', '3'), '0.1250')
        self.assertEqual(_av_weight('4', '4'), '0.0000')
        self.assertEqual(_av_weight('4', '9'), '0.0000')

    def test_readSTATION0(self):
        """Simple function to test the ability to read a test STATION0.HYP \
        file."""
        from eqcorrscan.utils.catalogue2DD import readSTATION0
        import os
        station_input_list = ['COVA', 'SOLU', 'drc', 'RPZ', 'NZ01', 'WHAT2',
                              'CAMEL', 'NZ20']
        STATION0_path = os.path.join(os.path.abspath(os.path.
                                                     dirname(__file__)),
                                     'test_data')
        station_output_list = readSTATION0(STATION0_path, station_input_list)
        # Check that the output file exists, and remove it
        self.assertTrue(os.path.isfile('station.dat'))
        for station_information in station_output_list:
            if station_information[0] == 'COVA':
                self.assertEqual(station_information,
                                 ('COVA', -43.6132, 169.968, 1477.0))
            elif station_information[0] == 'SOLU':
                self.assertEqual(station_information,
                                 ('SOLU', -43.90805, 169.60608333333334,
                                  1146.0))
            elif station_information[0] == 'drc':
                self.assertEqual(station_information,
                                 ('drc', -43.2411, 170.4515, 128.0))
            elif station_information[0] == 'RPZ':
                self.assertEqual(station_information,
                                 ('RPZ', -43.7164, 171.05388333333335, 453.0))
            elif station_information[0] == 'NZ01':
                self.assertEqual(station_information,
                                 ('NZ01', -44.498, 165.00198333333333,
                                  -4682.0))
            elif station_information[0] == 'WHAT2':
                self.assertEqual(station_information,
                                 ('WHAT2', -43.2793, 170.36038333333335, 95.0))
            elif station_information[0] == 'CAMEL':
                self.assertEqual(station_information,
                                 ('CAMEL', -42.874183333333335, 170.9622,
                                  60.0))
            elif station_information[0] == 'NZ20':
                self.assertEqual(station_information,
                                 ('NZ20', -40.2, 170.75, -714.0))
        # Check that things have been written correctly
        output_check_file = open('station.dat', 'r')
        for line in output_check_file:
            station = line[0:5].strip()
            latitude = float(line[5:].split()[0])
            longitude = float(line[5:].split()[1])
            depth = float(line[5:].split()[2].rstrip())
            station_information_out = (station, latitude, longitude, depth)
            # Check that this matches what we expect
            for station_information in station_output_list:
                if station_information[0] == station_information_out[0]:
                    self.assertEqual(round(station_information[1], 4),
                                     round(station_information_out[1], 4))
                    self.assertEqual(round(station_information[2], 4),
                                     round(station_information_out[2], 4))
                    self.assertEqual(round(station_information[3], 4),
                                     round(station_information_out[3] * 1000,
                                           4))
        output_check_file.close()
        os.remove('station.dat')

    def test_write_event(self):
        """Simple test function to test the writing of events.
        """
        from eqcorrscan.utils.catalogue2DD import sfiles_to_event
        from eqcorrscan.utils import Sfile_util
        from obspy import UTCDateTime
        import os
        import glob

        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_')
        sfile_list = glob.glob(os.path.join(testing_path, '*L.S??????'))
        event_list = sfiles_to_event(sfile_list)
        # Check that we have written a file
        self.assertTrue(os.path.isfile('event.dat'))
        f = open('event.dat', 'r')
        for line, event in zip(f, event_list):
            header = Sfile_util.readheader(event[1])
            event_id_input = event[0]
            output_event_info = line.strip().split()
            # Check that the event id's match
            self.assertEqual(event_id_input, int(output_event_info[-1]))
            time_string = str(header.origins[0].time.year) +\
                str(header.origins[0].time.month).zfill(2) +\
                str(header.origins[0].time.day).zfill(2)+'  ' +\
                str(header.origins[0].time.hour).rjust(2) +\
                str(header.origins[0].time.minute).zfill(2) +\
                str(header.origins[0].time.second).zfill(2) +\
                str(header.origins[0].time.microsecond)[0:2].zfill(2)
            self.assertEqual(output_event_info[0:2], time_string.split())
            self.assertEqual(header.origins[0].latitude,
                             float(output_event_info[2]))
            self.assertEqual(header.origins[0].longitude,
                             float(output_event_info[3]))
            self.assertEqual(header.origins[0].depth / 1000,
                             float(output_event_info[4]))
            if header.magnitudes[0]:
                self.assertEqual(header.magnitudes[0].mag,
                                 float(output_event_info[5]))
            if header.origins[0].time_errors.Time_Residual_RMS:
                self.assertEqual(header.origins[0].time_errors.
                                 Time_Residual_RMS,
                                 float(output_event_info[-2]))
        os.remove('event.dat')

if __name__ == '__main__':
    unittest.main()
