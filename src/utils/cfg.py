def read_cfg(cfg_path):
    import configparser

    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)

    return cfg


def get_cfg_items(cfg, section_name):
    # 获取所有的配置项
    items = cfg.items(section_name)
    return dict(items)


if __name__ == '__main__':
    config = read_cfg('data/cfg.ini')

    items_dict = get_cfg_items(config, 'digest')
    print(items_dict)

    # # 读取cfg文件
    # config = read_cfg('data/cfg.ini')
    #
    # paths_dict = get_cfg_items(config, "paths")

