#define L1_DATA_PACKET 1 /* Layer 1 directive for data/command packets */ 


#define L1_ACK_PACKET 2 /* Layer 1 directive for ACK packets */ 


#define L1_NAK_PACKET 3 /* Layer 1 directive for NAK packets */ 




#define SIZE_OF_HOST_POWER_SPECTRUM 60 /* 60 values for 0.5 hz - 30.0 hz */


#define SPECTRA_NUMOFCHAN 2 /* Number of spectra per message */


#define SIZE_M_SPECTRA ((SIZE_OF_HOST_POWER_SPECTRUM + 1)*4)


#define SIZE_M_PROCESSED_VARS 120 


#define SIZE_M_PROCESSED_VARS_AND_SPECTRA (SIZE_M_PROCESSED_VARS + SIZE_M_SPECTRA)


#define SIZE_M_PROCESSED_VARS_WITH_EXTRA_VARS 156 


#define SIZE_M_PROCESSED_VARS_AND_SPECTRA_WITH_EXTRA_VARS (SIZE_M_PROCESSED_VARS_WITH_EXTRA_VARS + SIZE_M_SPECTRA)


#define SIZE_M_PROCESSED_VARS_LABELS 120


#define SIZE_M_PROCESSED_EXTRA_VARS_LABELS 192


#define PERCENT_DIVISOR 10
#define DB_POWER_DIVISOR 100
#define FREQUENCY_DIVISOR 100


#define IMPED_MSG_DATA_LEN 49 /* Lengths of message-dependent data in bytes */


#define ERROR_MSG_DATA_LEN 62
#define REVISION_INFO_DATA_LEN 28
#define EEG_SNIPPET_HEADER_DATA_LEN 28


#define EEG_SNIPPET_PROCESSED_DATA_LEN 40
#define EEG_SNIPPET_PROCESSED_DATA_WITH_EXTRA_VARS_LEN 64


#define BIS_HISTORY_HEADER_DATA_LEN 24


#define MAX_BIS_HISTORY_DATA_LEN 1020


#define MAX_BIS_HISTORY_DATA_RECORDS 85 /* Max records per big History message */


#define EVENT_MSG_DATA_LEN 31 


#define START_CASE_RECORD 1
#define DATA_RECORD 2
#define EVENT_RECORD 3


#define CLOCK_ADJUST_RECORD 5


#define TIME_RECORD 6 

struct processed_vars_and_spectra_msg

{ struct processed_vars_msg processed_vars;
struct spectra_msg spectra;


};
 

struct processed_vars_and_spectra_with_extra_vars_msg

{ struct processed_vars_with_extra_vars_msg processed_vars_extra;
struct spectra_msg spectra;


};
 

struct processed_vars_msg{ struct dsc_info_struct proc_dsc_info;
 /* Info. on dsc, pic and selftest */ 

/* Impedance info for 2 channels. */

struct impedance_info_struct impedance_info[2];
 



/* Host filter settings:

Byte 0 (LSByte) - High pass filter setting

0 - 0.25 Hz 

1 - 1.00 Hz 

2 - 2.00 Hz 



3 - 2.5 Hz 

Byte 1 - Low pass filter setting

0 - none 

1 - 30 Hz 

2 - 50 Hz 

3 - 70 Hz 



Byte 2 - Notch filter setting

0 - none 

1 - 50 Hz 

2 - 60 Hz 

3 - 50 & 60 Hz */



unsigned long host_filt_setting;
 

/* Spectral and Bispectral smoothing rates:

Byte 0 (LSByte) - Spectral smoothing rate

0 ? 0 seconds 

1 ? 5 seconds 

2 ? 10 seconds (default)

3 ? 30 seconds 

4 ? 60 seconds 



Byte 1 - Bispectral smoothing rate

0 ? 15 seconds (default)

1 ? 30 seconds 

2 ? 60 seconds 

3 ? 10 seconds */



unsigned long host_smoothing_setting;
 

/* Spectral and Bispectral artifact detection masks. */

unsigned long host_spectral_art_mask;


unsigned long host_bispectral_art_mask;
 



/*--------------------------------------------------------------+

| trend variables for 3 channels (ch1, ch2 and ch12 combined) |

+--------------------------------------------------------------*/ 



struct be_trend_variables_info trend_variables[3];
};
 

struct processed_vars_with_extra_vars_msg{ struct dsc_info_struct proc_dsc_info;
 /* Info. on dsc, pic and selftest */ 

/* Impedance info for 2 channels. */

struct impedance_info_struct impedance_info[2];
 



/* Host filter settings:

Byte 0 (LSByte) - High pass filter setting

0 - 0.25 Hz 

1 - 1.00 Hz 

2 - 2.00 Hz 

3 - 2.5 Hz 



Byte 1 - Low pass filter setting

0 - none 

1 - 30 Hz 

2 - 50 Hz 

3 - 70 Hz 



Byte 2 - Notch filter setting

0 - none 

1 - 50 Hz 

2 - 60 Hz 

3 - 50 & 60 Hz */



unsigned long host_filt_setting;
 

/* Spectral and Bispectral smoothing rates:

Byte 0 (LSByte) - Spectral smoothing rate 



0 ? 0 seconds 1 ? 5 seconds 2 ? 10 seconds (default)3 ? 30 seconds 4 ? 60 seconds 

Byte 1 - Bispectral smoothing rate

0 ? 15 seconds (default)

1 ? 30 seconds 

2 ? 60 seconds 

3 ? 10 seconds */



unsigned long host_smoothing_setting;
 

/* Spectral and Bispectral artifact detection masks. */

unsigned long host_spectral_art_mask;


unsigned long host_bispectral_art_mask;
 



/*--------------------------------------------------------------+

| trend variables for 3 channels (ch1, ch2 and ch12 combined) |

| with extra variables |

+--------------------------------------------------------------*/ 



struct be_trend_variables_extra_info trend_variables_extra[3];
};
 

struct dsc_info_struct

{ unsigned char dsc_id;
 /* DSC ID from status nibble 1 */unsigned char dsc_id_legal;
 /* Flag when non-zero indicates that a legal

dsc is connected */

unsigned char pic_id;
 /* (Sensor Type * 10)+ PIC ID */ 




#define SINGLE_CHANNEL_OR_SENSOR_TYPE 1 /* Sensor Plus */
#define DUAL_CHANNEL_OR_SENSOR_TYPE 2 /* Quatro */
#define PEDIATRIC_OR_SENSOR_TYPE 3 /* Pediatric */
#define DUAL_CHANNEL_ICU_SENSOR_TYPE 4 /* Extend */
#define DEMO_SENSOR_TYPE 5 /* Demo (Plus) */
#define PEDIATRIC_XP_SENSOR_TYPE 7 /* Pediatric XP */
#define SINGLE_CHANNEL_SENSOR_SIMULATOR 8 /* Plus Simulator */
#define DUAL_CHANNEL_SENSOR_SIMULATOR 9 /* Quatro Sim. */
#define SEMI_REUSABLE_SENSOR_TYPE    12 /* SRS Type 1*/
#define EXTEND_DEMO_SENSOR_TYPE 13 /* Demo (Extend) */
#define BILATERAL_SENSOR_TYPE 14 /* Bilateral */
#define BILATERAL_SENSOR_TYPE_1 16 /* new Bilateral */
#define FOUR_CHANNEL_SENSOR_SIMULATOR 17 /* 4 channel Sim. */ 

unsigned char pic_id_legal;
 /* if non-zero, a legal pic is connected. */unsigned short dsc_numofchan;
 /* Number of channels on the DSC connected.Valid only when dsc_id is legal. */unsigned short quick_test_result;
 /* If zero, test passed.If non-zero, the bit fields in theresult indicate which test(s) failed. */ 


#define QUICK_SELFTEST_PASS 0 


#define QUICK_GAIN_TEST_BIT 0x1 /* Set when avg noise test fails. */


#define QUICK_NOISE_TEST_BIT 0x2 /* Set when blocked gain test fails. */


#define QUICK_TEST_FAIL_BIT 0x4 /* Set for timeout, DSC disconnect and



other failures during the test .*/


/* DSC gain (in uV/ADC) = dsc_gain_num/dsc_gain_divisor */

signed long dsc_gain_num;


signed long dsc_gain_divisor;
 



/* DSC offset (in uV/ADC) = dsc_offset_num/dsc_offset_divisor */

signed long dsc_offset_num;


signed long dsc_offset_divisor;




};
 

struct impedance_info_struct

{ /*----------------------------------------------------------+| Impedance in units of 100 Ohms (or 1000 Ohms for ground || impedances). |+----------------------------------------------------------*/


#define GND_SCALE 1000 
#define NON_GND_SCALE 100 

/* Limits for the impedance values */
#define MAX_GND_VALUE 100.0 /* in Kohms */
#define MAX_COM_VALUE 150.0 /* in 100 ohms */
#define MAX_POS_VALUE 75.0 /* in 100 ohms */ 

unsigned short impedance_value;
 /* This is the value of the

impedance for the channel */



unsigned short imped_test_result;
 /* Result of the impedance test.

If zero, test passed.

If non-zero, the bit fields in the



result indicate which test(s) failed. */ 


#define IMPED_TEST_PASS 0 /* Test passed */


#define IMPED_TEST_FAIL 0x1 /* High impedance */


#define IMPED_TEST_FAIL_CLIP 0x2 /* Noise */


#define IMPED_TEST_FAIL_LEADOFF 0x4 /* Lead off */

};
 



struct be_trend_variables_info

{ 



signed short burst_suppress_ratio;
 /* index variable giving percent ofsuppressed seconds in last 63 sec.for selected channel. range from0 - 1000 in .1% steps */

signed short spectral_edge_95;
 /* in HZ ranged from 0-30.0 Hz in

units of 0.01 Hz */

signed short bis_bits;
 /* BIS field debug data */

signed short bispectral_index;
 /* Ranges from 0 - 100 */

signed short bispectral_alternate_index;
 /* same as above */

signed short bispectral_alternate2_index;
 /* same as above */

signed short total_power;
 /* in dB with respect to .01 uV rms. ranged



from 0 to 100 dB in 0.01 units */

signed short emg_low;
 /* in dB with respect to .01 uV rms. ranged

from 0 to 100 dB in 0.01 units */



signed long bis_signal_quality;
 /* index variable giving the signal

quality of the bisIndex which is

combined with BSR 0 - 1000 in .1% 

steps */



unsigned long second_artifact;
 /* bit field indicating TYPE OF ARTIFACTfor the last second. */};
 

struct be_trend_variables_extra_info{ 

signed short burst_suppress_ratio;
 /* index variable giving percent ofsuppressed seconds in last 63 sec.for selected channel. range from 

0 - 1000 in .1% steps */

signed short spectral_edge_95;
 /* in HZ ranged from 0-30.0 Hz in



units of 0.01 Hz */

signed short bis_bits;
 /* BIS field debug data */

signed short bispectral_index;
 /* Ranges from 0 - 100 */

signed short bispectral_alternate_index;
 /* same as above */

signed short bispectral_alternate2_index;
 /* same as above */

signed short total_power;
 /* in dB with respect to .01 uV rms. ranged



from 0 to 100 dB in 0.01 units */

signed short emg_low;
 /* in dB with respect to .01 uV rms. ranged

from 0 to 100 dB in 0.01 units */



signed long bis_signal_quality;
 /* index variable giving the signal

quality of the bisIndex which is

combined with BSR 0 - 1000 in .1% 

steps */



unsigned long second_artifact;
 /* bit field indicating TYPE OF ARTIFACT

for the last second. */

signed short burst_per_min;
 /* ranges from 0 to 30 */

signed short rfu1;
 /* reserved for future use */

signed short rfu2;


signed short rfu3;


signed short rfu4;


signed short rfu5;




};
 

struct spectra_msg

{ unsigned short spect_numofchan;
unsigned short spect_size;
signed short power_spectrum[SPECTRA_NUMOFCHAN][SIZE_OF_HOST_POWER_SPECTRUM];


};
 

struct revision_info_msg

{ char system_revision[4];
 /* Binary bytes: first is major rev. number, */char host_revision[4];
 /* second is minor rev. number, third is */char bis_eng_revision[4];
 /* not defined, fourth is always 0. */char protocol_revision[4];
 /* Displayed as major.minor */char boot_revision[4];
char hardware_revision[4];
char serial_number[4];
 /* First 3 bytes make up a 24-bit unsigned */

};
 /* integer (LSB first), last byte is an alphabetic */

/* char (typ. ?A?). */

/* Displayed as char followed by integer. */ 



struct be_trend_variables_labels

{ char burst_suppress_ratio_label[12];
 /* ASCII strings */char spectral_edge_95_label[12];
char bis_bits_label[12];
char bispectral_index_label[12];
char bispectral_alternate_index_label[12];
char bispectral_alternate2_index_label[12];
char total_power_label[12];
char emg_low_label[12];
char bis_signal_quality_label[12];
char second_artifact_label[12];


};
 

struct be_trend_variables_extra_labels{ char burst_suppress_ratio_label[12];
 /* ASCII strings */ 

char spectral_edge_95_label[12];


char bis_bits_label[12];


char bispectral_index_label[12];


char bispectral_alternate_index_label[12];


char bispectral_alternate2_index_label[12];


char total_power_label[12];


char emg_low_label[12];


char bis_signal_quality_label[12];


char second_artifact_label[12];


char burst_per_min_label[12];


char rfu1_label[12];


char rfu2_label[12];


char rfu3_label[12];


char rfu4_label[12];


char rfu5_label[12];




};
 

struct time_date

{ short second;
 /* second, 0 - 59 */short minute;
 /* minute, 0 - 59 */short hour;
 /* hour, 0 - 23 */short day;
 /* day, 1 - 31 */short month;
 /* month, 1 - 12 */short year;
 /* 4-digit year (with century) */

};
 

struct snippet_info_msg

{ struct time_date event_date_time;
 /* Date/Time */char system_revision[4];
 /* System revision number */char serial_number[4];
 /* Monitor serial number */short num_raw_records;
 /* Number of raw data records */short num_proc_records;
 /* Number of processed data records */unsigned char dsc_id;
 /* DSC ID */unsigned char pic_id;
 /* PIC ID */short padding;
 /* fills out last word */

};
 

struct snippet_processed_vars_msg{ struct snippet_trend_variables_info snippet_vars[2];
};
 

struct snippet_trend_variables_info

{ short impedance_value;
short burst_suppress_ratio;
short bis_bits;
short bispectral_index;
short bispectral_alternate_index;
short bispectral_alternate2_index;
short emg_low;
short bis_signal_quality;
unsigned long second_artifact;


};
 

struct snippet_processed_extra_vars_msg{ struct snippet_trend_extra_variables_info snippet_vars[2];
};
 

struct snippet_trend_extra_variables_info

{ short impedance_value;
short burst_suppress_ratio;
short bis_bits;
short bispectral_index;
short bispectral_alternate_index;
short bispectral_alternate2_index;
short emg_low;
short bis_signal_quality;
unsigned long second_artifact;
short burst_per_min;


short rfu_1;
short rfu_2;
short rfu_3;
short rfu_4;
short rfu_5;
};
  

struct history_info_msg{ struct time_date date_time;
char system_revision[4];
char serial_number[4];
long num_records;
};
  /* Date/Time *//* System revision number *//* Monitor serial number *//* Number of data records */  

struct history_data_msg{ unsigned long num_records;
 



struct bh_record_struct data[MAX_BIS_HISTORY_DATA_RECORDS];


};
 

struct bh_record_struct{ 

union bh_record_body {struct bh_startcase_or_clockadjust start;
struct bh_data data;


} body;
};
 

struct bh_startcase_or_clockadjust

{ unsigned char type;
 /* Record type (start case = 1, clock adjust = 5) */unsigned char pad1;
 /* padding */short pad2;
 /* padding */long time;
 /* Time of start case or new clock time, as utime */long case_id;
 /* Case ID */long pad3;
   /* padding */ 

};
 

struct bh_data

{ unsigned char type;
 /* Record type (data = 2) */unsigned char avg_bis;
 /* Average BIS for 1 minute */unsigned char min_bis;
 /* Minimum BIS for 1 minute */unsigned char max_bis;
 /* Maximum BIS for 1 minute */unsigned char avg_bisalt;
 /* Average BISALT for 1 minute */unsigned char avg_bisalt2;
 /* Average BISALT2 for 1 minute */unsigned char avg_sqi;
 /* Average SQI for 1 minute */ 

unsigned char avg_emg;
 /* Average EMG for 1 minute */

unsigned char avg_sr;
 /* Average SR for 1 minute */

unsigned char avg_imped[2];
 /* Average combined impedances for 1 minute */

short bis_bits;
 /* BIS Bits from last update in minute */

byte avg_burst_per_min;
 /* Unused if extra vars have not been requested */

byte reserved;
 /* Reserved for future use */



};
 

/* VISTA Binary Mode Structures */ 


#define SIZE_M_PROCESSED_VARS_4B 288 


#define SIZE_M_PROCESSED_VARS_AND_SPECTRA_4B (SIZE_M_PROCESSED_VARS_4B + SIZE_M_SPECTRA) 



struct revision_info_msg 



{ 	/* Binary bytes: first is major rev. number, */ 

/* second is minor rev. number 

/* third and fourth are always 0. */ 

/* Displayed as major.minor */ 



char system_revision[4];
 /* VISTA Master Software Revision */ 

char host_revision[4];
 /* VISTA Application Revision */ 

char bis_eng_revision[4];
 /* BISx Software Revision */ 

char protocol_revision[4];
 /* Serial Protocol Revision */ 

char boot_revision[4];
 /* BISx Serial Number ? see serial_number below 



BISx last digit is ?B? 



BISx4 last digit is ?4? */ 

char hardware_revision[4];
 // VISTA Hardware Revision 

char serial_number[4];
 /* First 3 bytes make up a 24-bit unsigned */ 



/* integer (LSB first), last byte is first alphabetic */ 

/* char. Displayed as char followed by integer. */ 

/* e.g. VT01234 -> 0xD2 0x04 0x00 0x56 */ 

/* VISTA: last byte is ?V? */ 



};
 



struct be_bilateral_trend_variables_labels

{ char burst_suppress_ratio_label[12];
char spectral_edge_95_label[12];
char spectral_edge_50_label[12];
char bis_bits_label[12];
char bispectral_index_label[12];
char bispectral_alternate_index_label[12];
char bispectral_alternate2_index_label[12];
char total_power_label[12];
char emg_low_label[12];
char bis_signal_quality_label[12];
char second_artifact_label[12];
char burst_per_min_label[6];
char reserved_label[6];
char asym_label[12];
char std_bis_label[12];
char std_emg_label[12];
char reserved_label_0[12];
char reserved_label_1[12];
char reserved_label_2[12];


};
 

struct processed_vars_and_spectra_msg_4b

{ struct processed_vars_msg_4b processed_vars;
struct spectra_msg spectra;


};
 

struct processed_vars_msg_4b{ struct dsc_bilateral_info_struct proc_dsc_info;
 /* Info. on dsc, pic and selftest */ 

/* Impedance info for up to 4 channels. */

struct impedance_info_struct impedance_info[4];
 



/* Host filter settings:

Byte 0 (LSByte) - High pass filter setting

0 - 0.25 Hz 

1 - 1.00 Hz 

2 - 2.00 Hz 

3 - 2.5 Hz 



Byte 1 - Low pass filter setting

0 - none 

1 - 30 Hz 

2 - 50 Hz 

3 - 70 Hz 



Byte 2 - Notch filter setting

0 - none 

1 - 50 Hz 

2 - 60 Hz 

3 - 50 & 60 Hz 



Byte 3 - Operating environment: OPERATING_ENV_OR or OPERATING_ENV_ICU

*/

unsigned long host_filt_setting;
 



/* Spectral and Bispectral smoothing rates:

Byte 0 (LSByte) - Spectral smoothing rate

Byte 1 - Bispectral smoothing rate



*/

unsigned long host_smoothing_setting;


/* 



Spectral smoothing rates Bispectral smoothing rates

Option Rate (in sec) Option Rate (in sec)

0 0 0 15 

1 5 1 30 

2 10 2 60 

3 30 3 10 

4 60 

*/ 




# define MIN_SPECT_SMOOTHING_RATE 0 
# define MAX_SPECT_SMOOTHING_RATE 4 
# define MIN_BISPECT_SMOOTHING_RATE 0 
# define MAX_BISPECT_SMOOTHING_RATE 3 

/* Spectral and Bispectral artifact detection masks.

These are provided for validation purposes only and should be ignored

during normal monitoring. */ 



unsigned long host_spectral_art_mask;


unsigned long host_bispectral_art_mask;
 



/*------------------------------------------------------------------+

| trend variables for up to 4 channels (ch1, ch2, ch3, ch4)

+------------------------------------------------------------------*/ 



struct be_bilateral_trend_variables_info trend_variables[4];
 

unsigned short sqi_left_index;
 /* Index into trend_variables for selecting

the 'bis_signal_quality' variable to use for

sqi threshold comparisons.

Possible values: 0 - 4 */



unsigned short sqi_right_index;
 /* Index into trend_variables for selecting

the 'bis_signal_quality' variable to use for

sqi threshold comparisons.

Possible values: 0 - 4 */



unsigned short bis_left_index;
 /* Index into trend_variables for selectingthe 'bispectral_index' variable to display as BIS.Possible values: 0 - 4 */

unsigned short bis_right_index;
 /* Index into trend_variables for selecting 

the 'bispectral_index' variable to display as BIS.Possible values: 0 - 4 */};
 

struct dsc_bilateral_info_struct{ unsigned char dsc_id;
 /* DSC ID from status nibble 1 */ 


#define DSC_NONE_ID 0 /* No DSC connected */
#define DSC_BISX4_ID 13 

/* Only DSC IDs of DSC_NONE_ID or DSC_BISX4_ID are considered legal;


For all other values, dsc_id_legal is set to zero to indicate

illegal DSC ID. */ 



unsigned char dsc_id_legal;
 /* Flag when non-zero indicates that a legal

dsc is connected */

unsigned char pic_id;
 /* Sensor/pigtail ID */ 




#define PIC_NONE_ID 0 
#define PIC_SS_ID 7 

/* Legal PIC IDs are:

PIC_NONE_ID

PIC_SS_ID + (Sensor Type * 10)

For all other IDs, pic_id_legal is set to zero. */ 



unsigned char pic_id_legal;
 /* if non-zero, a legal pic is connected. */

unsigned short dsc_numofchan;
 /* Number of channels on the DSC connected.

Valid only when dsc_id is legal. */

unsigned short quick_test_result;
 /* If zero, test passed.

If non-zero, the bit fields in the

result indicate which test(s) failed. */ 



other failures during the test .*/


#define QUICK_TEST_RESULT_VALID_BIT 0x8 /* Set when the quick test has notbeen run at all. Cleared when the test is done atleast once. */ 

unsigned short dsc_update_status;
 /* Bit field indicating DSC update status. */
#define DSC_UPDATE_DONE 0 
#define DSC_UPDATE_IN_PROGRESS 0x1 
#define DSC_UPDATE_FAILED 0x2 

unsigned short pad1;
 // for alignment 

/* DSC gain (in uV/ADC) = dsc_gain_num/dsc_gain_divisor */

signed long dsc_gain_num;


signed long dsc_gain_divisor;
 



/* DSC offset (in uV/ADC) = dsc_offset_num/dsc_offset_divisor */

signed long dsc_offset_num;


signed long dsc_offset_divisor;
 



/****************/

/* Sensor Info */

/****************/ 



unsigned short sensor_status;
 /* high nibble is status nibble 5 */ 


#define SENSOR_VALID 0 
#define SENSOR_INVALID 1 
#define TOO_MANY_USES 2 
#define SENSOR_EXPIRED 3 

/* Nibble-5 status Bit 
# */
#define SENSOR_POWERED_OFF 12 
#define SENSOR_OVERCURRENT_BIT 13 
#define SENSOR_POS_GND_FAULT_BIT 14 
#define SENSOR_NEG_GND_FAULT_BIT 15 


#define SENSOR_UNSUPPORTED 0x20 /* Bit 5: Unsupported sensor */


#define SRS_ELECTRODES_NOT_CONNECTED 0x40 /* Bit 6: 0 if electrodesdetected, 1 if notdetected. */


#define NEW_VALIDITY_UNKNOWN 0x3f 
#define EXTEND_TYPE_SENSOR 0x200 /* Bit 9: Sensor activates

Burst Count features */
#define SENSOR_SIMULATOR 0x400 /* Bit 10: BIS values are blanked */
#define DEMO_DEVICE 0x800 /* Bit 11: Must display Demo

Device message */ 


#define NEW_SS_BITS_MASK 	(0xffff080 | SRS_ELECTRODES_NOT_CONNECTED |\SRS_TYPE_SENSOR |\EXTEND_TYPE_SENSOR |\SENSOR_SIMULATOR |\DEMO_DEVICE)

struct bilateral_sensor_description_struct sensor_desc;


unsigned short pad2;
 // for alignment

unsigned long lot_code;


unsigned short shelf_life;


unsigned short serial_number;


unsigned long usage_count;




};
 

struct bilateral_sensor_description_struct

{ unsigned char sensor_name[12];
/* 12 character alpabetic string uniquely identifying the sensor. */ 

unsigned char sensor_type;
 /* Unique sensor identifier. */ 

unsigned char sensor_graphic_type;


/* Number of channels to display during impedance checking.

If sensor_graphic_type is 6, impedance values should be checked for

channel 1,channel 2, channel 3 and channel 4.

*/ 



unsigned char sensor_eeg_channels;


/* Number of channels of filtered EEG to be displayed. */ 



unsigned char eeg_left_display_index;


unsigned char eeg_right_display_index;


/* This index indicates the channel to be displayed. */ 



unsigned char sensor_sqi_channels;


/* Number of channels of SQI values to be displayed. */ 



unsigned char sqi_left_display_index;


unsigned char sqi_right_display_index;


/* Use this as an index into the



channel data structure with sqi value. */ 

unsigned char sensor_emg_channels;
 /* Number of channels of EMG values to be displayed. */ 

unsigned char emg_left_display_index;


unsigned char emg_right_display_index;
 



unsigned char pad1;


/* Use this as an index into the

channel data structure with emg value. */ 



unsigned char case_id[4];


/* Alphanumeric characters representing case ID formed by combining

sensor serial number and BIS Engine serial number.

Display format: xxxx



*/};
 

struct be_bilateral_trend_variables_info{ 

signed short burst_suppress_ratio;
 /* index variable giving percent of

     suppressed seconds in last 63 sec. 

     for selected channel. range from 

0 - 1000 in .1% steps */



signed short spectral_edge_95;
 /* in HZ ranged from 0-30.0 Hz in

units of 0.01 Hz */

signed short spectral_edge_50;
 /* in HZ ranged from 0-30.0 Hz in



units of 0.01 Hz */

signed short bis_bits;
 /* BIS field debug data */

signed short bispectral_index;
 /* Ranges from 0 - 100 */

signed short bispectral_alternate_index;
 /* same as above */

signed short bispectral_alternate2_index;
 /* same as above */

signed short total_power;
 /* in dB with respect to .01 uv rms. ranged



from 0 to 100 dB in 0.01 units */

signed short emg_low;
 /* in dB with respect to .01 uv rms. ranged



from 0 to 100 dB in 0.01 units */

unsigned short pad1;
 /* for alignment */

signed long bis_signal_quality;
 /* index variable giving the signal



quality of the bisIndex which iscombined with BSR 0 - 1000 in .1% steps */

unsigned long second_artifact;
 /* bit field indicating TYPE OF ARTIFACT

for the last second. */ 



unsigned char burst_count;
 /* bursts/minute 0 - 30 */
#define BURST_COUNT_MASK 0x3f /* mask for burst count value */
#define BURST_COUNT_BLANK 0x40 /* burst count blanked if this bit is set */
#define BURST_COUNT_ENABLE 0x80 /* burst count enabled if this bit is set */

unsigned char reserved_byte_var;
 /* reserved*/

signed short asym_index;


signed short std_bis;
 /* standard deviation of BIS */

signed short std_emg;
 /* standard deviation of EMG */

signed short reserved_short_var_0;


signed short reserved_short_var_1;


signed short reserved_short_var_2;


unsigned short pad2;
 /* for alignment */



};
