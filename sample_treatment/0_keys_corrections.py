# import re
#
# def replace_dict_keys_everywhere(input_file, output_file, replacement_dict):
#     # Sort keys by length descending to prevent partial replacements (e.g., O3_500 before O3_5008A)
#     sorted_keys = sorted(replacement_dict.keys(), key=len, reverse=True)
#
#     # Build pattern using alternation, escaping keys properly
#     pattern = re.compile('|'.join(map(re.escape, sorted_keys)))
#
#     def replacer(match):
#         return replacement_dict[match.group(0)]
#
#     with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
#         for line in f_in:
#             new_line = pattern.sub(replacer, line)
#             f_out.write(new_line)
#
#
#
#
# conv_dict = {'O2_3727A': 'O2_3726A', 'O2_3730A': 'O2_3729A', 'Ne3_3870A': 'Ne3_3869A',
#              'H1_3890A': 'H1_3889A', 'He1_3890A': 'He1_3889A', 'H1_4103A': 'H1_4102A',
#              'H1_4342A': 'H1_4340A', 'O3_4364A': 'O3_4363A', 'H1_4863A': 'H1_4861A',
#              'O3_4960A': 'O3_4959A', 'O3_5008A': 'O3_5007A', 'He1_5877A': 'He1_5876A',
#              'N2_6550A': 'N2_6548A', 'H1_6565A': 'H1_6563A', 'N2_6585A': 'N2_6583A',
#              'S2_6718A': 'S2_6716A', 'S2_6733A': 'S2_6731A', 'S3_9071A': 'S3_9068A',
#              'S3_9533A': 'S3_9530A'}
#
# # Load configuration
# cfg_file = '../CAPERS_v3.toml'
# new_file = "../new_version.toml"
# replace_dict_keys_everywhere(cfg_file, new_file, conv_dict)


import lime
# Load configuration
cfg_file = '../CAPERS_v3.toml'
capers_cfg = lime.load_cfg(cfg_file)

lines_list = capers_cfg['lines_data']['lines_database']
bands = lime.line_bands(line_list=lines_list, vacuum_waves=True)

lime.save_frame('/home/vital/Dropbox/Astrophysics/Data/CAPERS/source/CAPERs_lines_database.txt', bands)