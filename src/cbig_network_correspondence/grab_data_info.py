"""Module providing argument parser"""
import configparser


class MyClass:
    def __init__(self):
        self.category = None
        self.name = None
        self.space = None
        self.type = None
        self.netassign = None
        self.threshold = None


def str2bool(str_val):
    """
    Convert string to boolean
    """
    if str_val is None:
        return False
    return str_val.lower() in ("yes", "true", "t", "1")


def str2num(str_val):
    """
    Convert string to number
    """
    if str_val is None:
        return None
    if str_val.lower() == 'inf':
        return float('inf')
    if str_val.lower() == '-inf':
        return float('-inf')
    try:
        return int(str_val)
    except ValueError:
        return float(str_val)


def thresh_str_to_num(str_val):
    """
    Convert threshold range from string to number
    """
    if str_val is None:
        return None

    th_l, th_h = map(str2num, str_val[1:-1].replace(',', ' ').split())
    return th_l, th_h


def config_get(config, variable_name):
    """
    Reads in variable from the config. If variable does not exist, return None.
    """
    try:
        value = config.get("data_info", variable_name)
    except configparser.NoOptionError:
        value = None
    return value


def read_config(config_file):
    """
    Reads in data from a config file.

    :param config_file: path to the config file.
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    data_category = config_get(config, "Data_Category")
    data_name = config_get(config, "Data_Name")
    data_space = config_get(config, "Data_Space")
    data_type = config_get(config, "Data_Type")
    data_network_assignment = config_get(config, "Data_NetworkAssignment")
    data_threshold = config_get(config, "Data_Threshold")

    args = MyClass()
    args.category = data_category
    args.name = data_name
    args.space = data_space
    args.type = data_type
    args.netassign = str2bool(data_network_assignment)
    args.threshold = thresh_str_to_num(data_threshold)

    args = check_args_format(args)

    return args


def check_args_format(args):
    """
    Checks the format of the arguments.

    args: arguments.
    """
    if args.name is None:
        raise ValueError("Data_Name is not specified.")
    if args.space is None:
        raise ValueError("Data_Space is not specified.")
    if args.type is None:
        raise ValueError("Data_Type is not specified." +
                         "Please specify 'Hard' for hard parcellation." +
                         "Please specify 'Soft' for soft parcellaiton." +
                         "Please specify 'Metric' for metric data," +
                         "e.g., probability map, task contrast map.")
    if args.netassign is None:
        print("Data_NetworkAssignment is not specified.")
        args.netassign = False
    if args.threshold is None:
        if args.type == 'Soft':
            args.threshold = thresh_str_to_num('[0, Inf]')
            print("Data threshold is not specified for soft parcellation" +
                "or metric data. Will use the default threshold [0, Inf].")
    return args
