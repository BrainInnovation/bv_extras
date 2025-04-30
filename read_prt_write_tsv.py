"""
Save PRT as TSV file
last changed 2025-04-30

"""

import numpy as np
import bvbabel
from PyQt5.QtWidgets import QFileDialog, QApplication, QMessageBox, QInputDialog
import sys
import os
import traceback


def get_tr():
    """Prompt user for time to repeat (TR) in milliseconds."""
    while True:
        tr, ok = QInputDialog.getInt(None, "TR (ms)", "Enter the TR (100–10000):")
        if not ok:
            QMessageBox.information(None, "Canceled", "Operation canceled.")
            sys.exit()
        if 100 <= tr <= 10000:
            return tr
        QMessageBox.warning(None, "Invalid", "Enter a value between 100 and 10000.")


def convert_to_tsv(prt_path, tr):
    """Convert a PRT file to a BIDS-style TSV file."""
    print("Starting conversion...")
    print(f"Input file: {prt_path}")
    print(f"TR: {tr} ms")

    header, data = bvbabel.prt.read_prt(prt_path)
    rows = []

    for cond in data:
        start = np.array(cond['Time start'], dtype=int)
        stop = np.array(cond['Time stop'], dtype=int)

        if header['ResolutionOfTime'] != 'msec':
            start = (start - 1) * tr / 1000
            stop = stop * tr / 1000
        else:
            start = start / 1000
            stop = stop / 1000

        duration = np.round(stop - start, 4)
        condition = cond['NameOfCondition']
        rows.extend(zip(start, duration, [condition] * len(start)))

    output_path = f"{os.path.splitext(prt_path)[0]}_events.tsv"
    np.savetxt(output_path, rows, fmt='%s', delimiter='\t',
               header='onset\tduration\ttrial_type', comments='')

    print(f"Output file: {output_path}")
    print("Conversion successful.")
    QMessageBox.information(None, "Done", f"Saved TSV: {output_path}")


def main():
    app = QApplication(sys.argv)

    try:
        prt_path, _ = QFileDialog.getOpenFileName(None, "Select PRT File", os.path.expanduser("~"), "Protocol (*.prt)")
        if not prt_path:
            QMessageBox.warning(None, "No File", "No file selected.")
            sys.exit()

        tr = get_tr()
        convert_to_tsv(prt_path, tr)
        app.quit()

    except Exception:
        QMessageBox.critical(None, "Error", "Unexpected error occurred.")
        traceback.print_exc()
        app.quit()

    sys.exit()


if __name__ == "__main__":
    main()
