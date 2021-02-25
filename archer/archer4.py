from __future__ import print_function
from archer.archer4_mw  import archer4_mw
from archer.archer4_visir import archer4_visir
import archer.utilities.fdeckToolbox as fdtbx


def archer4(image, attrib, first_guess, sector_info=None, para_fix=True, display_filename=None):

    # archer4:
    #
    # This function channels the operations into one of the archer4 variants, which have
    # the same I/O requirements but different internal logic.
    #
    # AJW, CIMSS, Apr 2020. *** If you have questions unfortunately I will not answer them. ***

    """
    
    This function inputs three dictionaries and two flags, and outputs three dictionaries, 
    as follows:

    ========== Input ==========================================================

    image: 2D grids corresponding to the image
        image['lat_grid']:  Latitude grid. See notes.
        image['lon_grid']:  Longitude grid
        image['data_grid']: Image data (brightness temp or 0-255 brightness value)
        *image['azm_grid']: Azimuthal direction of beam (available from VIISR and ATMS Level 1B)
        *image['zen_grid']: Zenith angle of beam (available from same)
            *: Optional. Used for highest accuracy parallax calculation if desired.

    attrib: Image attributes
        attrib['sat']: Satellite source (NOAA-16, Meteosat-10, F18, Aqua, GOES-16, etc). See notes.
        attrib['sensor']: Sensor instrument (SSMI, Imager, etc). See notes.
        attrib['scan_type']: 'Conical', 'Crosstrack' or 'Geo'. Used for parallax calculations
        attrib['archer_channel_type']: Name of channel according to archer rules. See notes.
        +attrib['nadir_lon']: Nadir longitude (East +).
            +: Necessary for Geo data only. Used for parallax calculation.

    first_guess: First guess of the center fix
        *first_guess['source']: Source of first guess estimate (fx, bt, NHCanfx, etc). See notes.
        *first_guess['time']: Time of observation, seconds from epoch at 1 Jan 1970.
        first_guess['vmax']: Estimated Vmax of the TC.
        first_guess['lat']: Estimated first guess latitute of TC center.
        first_guess['lon']: Estimated first guess longitude of TC center.
            *: Optional. Not used yet, but probably will be brought in later.

    sector_info: Dictionary available from GeoIPS (pyresample). Currently this is only used here
        for information passed on to the f-deck output string
        sector_info['storm_basin'] : Two-character basin code
        sector_info['storm_num'] : Two-digit storm code

    para_fix: Flag to control whether to apply parallax fix to the imagery. True [or False].

    display_filename: Path and filename of diagnostic display image. No image if None.


    ========== Output =========================================================

    in_dict: Variables as used in ARCHER calculations, presented here for diagnostics.
        in_dict['sensor'] = Same as attrib['archer_channel_type']
        in_dict['lon_mx'] = Navigation of image used here (either parallax fixed or not)
        in_dict['lat_mx'] = "
        in_dict['bt_mx'] = Image data, either in BT or pseudo-BT
        in_dict['time'] = Time of the image, not used by ARCHER.
        in_dict['op_lon'] = Same as first_guess['lon']
        in_dict['op_lat'] = Same as first_guess['lat']
        in_dict['op_vmax'] = Same as first_guess['vmax']
        in_dict['ring_weight'] = Relative weight of ring score (dependent on sensor)

    out_dict: ARCHER output.
        out_dict['archer_channel_type'] = Same as attrib['archer_channel_type']
        out_dict['ring_radius_deg'] = Radius of resolved inner eyewall, units: deg
        out_dict['score_by_radius_arr'] = Ring score as a function of radius
        out_dict['confidence_score'] = Metric of contour spacing near center fix.
        out_dict['alpha_parameter'] = Parameter of gamma distribution (uncertainty)
        out_dict['radius50percCertDeg'] = Radius of 50% certainty area for center fix
        out_dict['radius95percCertDeg'] = Radius of 95% certainty area for center fix
        out_dict['eye_prob'] = Empirically determined probability of an eye (85GHz and IR only)
        out_dict['center_lon'] = ARCHER center fix longitude
        out_dict['center_lat'] = ARCHER center fix latitude
        out_dict['weak_center_lon'] = Last-resort center fix that violates some rules
        out_dict['weak_center_lat'] = "

    score_dict: Intermediate products of ARCHER, for diagnostics and dependent algorithms.
        score_dict['lon_grid1'] = Resampled regular grid coordinates
        score_dict['lat_grid1'] = Resampled regular grid coordinates
        score_dict['data_grid1'] = Resampled image
        score_dict['spiral_score_grid'] = Grid representing center of spiral shearing 
        score_dict['penalty_grid'] = Grid that exercises light penalty for straying from first guess
        score_dict['ring_score_grid'] = Grid representing center of possible eye
        score_dict['ring_radius_grid'] = Grid of candidates for ring radius
        score_dict['ring_score_grid_full'] = 3D grid (lat x lon x radius) used for ERC apps
        score_dict['radial_gradient_4d'] = 4D grid (lat x lon x 2d of gradient) for diagnostics
        score_dict['fraction_input'] = Fraction of domain covered by real image data


    ========== Notes ===========================================================

    1. Image navigation rules:
        a. You must use the original navigation of the image in order to calculate parallax.
        b. Archer can calcuate parallax from a subsection (line and/or element) of a geo image
            or conical scan. However, to calculate parallax of crosstrack scans, Archer must have 
            a subsection with full lines (all the elements of the line). 
        c. Conical scans are expected to default to the format that is found in their Level 1B
            files, which is rows = lines and cols = elements. 
        d. *** Cross-track scans are expected to default to the opposite format, as found in *their*
            Level 1B files, which is rows = elements and cols = lines. For example, a grid of 
            ATMS data is expected to have 96 rows always, VIIRS DNB data is expected to have 
            4064 rows always, etc.
        e. The one exception to all this is if you include image['azm_grid'] and image['zen_grid'] 
            as inputs. Then archer can calculate parallax on any projection or subsection.


    2. attrib['sat']: Name of the satellite. This is not *yet* used in Archer, but may be
        necessary in the future. Follow the GeoIPS rules for satellite naming.


    3. attrib['sensor']: One of the following:
        For conical microwave sensors: ['SSMI', 'SSMIS', 'TMI', 'GMI', 'AMSRE', 'AMSR2']
        For crosstrack microwave sensors: ['AMSUB', 'MHS', 'ATMS']
        For polar imagers: ['VIIRS'] (Add AVHRR?)
        For geo imagers: ['Imager']


    4. attrib['archer_channel_type']: One of these valid strings: ['37GHz', '89GHz', '183GHz',
        'IR', 'SWIR', 'Vis', 'DNB']
        Used for directing channel-specific logic in archer.


    5. first_guess['source']: One of these strings:
        'fx': Forecast
        'bt': Best track
        'manual': Manual estimate
        'archer': Previous archer estimate
        'NHCfxan': Splice of NHC-generated analysis and forecast track
        This variable is for nothing yet, but it may be needed in the future.


    6. For more context on these methods, please refer to:

            Wimmers, A. J., and C. S. Velden, 2016: Advancements in objective multisatellite 
            tropical cyclone center fixing. J. Appl. Meteor. Climatol., 55, 197–212.

    """


    if 'GHz' in attrib['archer_channel_type']:

        in_dict, out_dict, score_dict = archer4_mw(
            image, attrib, first_guess, para_fix=para_fix, display_filename=display_filename)

    else:

        in_dict, out_dict, score_dict = archer4_visir(
            image, attrib, first_guess, para_fix=para_fix, display_filename=display_filename)

    # In the future, add an option for archer4_scat


    # Produce an fdeck-formatted string
    fdeck_str = fdtbx.generate_string(attrib, in_dict, out_dict, sector_info=sector_info)
    print('fdeck output:')
    print(fdeck_str)
    out_dict['fdeck_string'] = fdeck_str


    return in_dict, out_dict, score_dict
