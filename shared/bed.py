

class TruthStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


class BedWriter(object):
    def __init__(self, bed_fn):
        self.bed_fn = bed_fn
        self.bed_writer = open(self.bed_fn, 'w')
    def close(self):
        try:
            self.bed_writer.close()
        except:
            pass

    def write_row(self, ctg_name, start_pos, end_pos, extra_infos="", zero_index=True):
        if not zero_index:
            start_pos -= 1
            end_pos -= 1
        start_pos = str(start_pos)
        end_pos = str(end_pos)
        bed_format = '\t'.join([ctg_name, start_pos, end_pos]) + '\n'

        self.bed_writer.write(bed_format)
