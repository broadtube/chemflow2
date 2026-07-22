"""ビルトイン装置。

新しい装置はこのフォルダに 1 ファイル追加し、Unit を継承して residuals() を書く。
ここに import 追加すれば ``from chemflow2 import <NewUnit>`` で使えるようになる。
"""

from chemflow2.units.gibbs import CanteraError, GibbsReactor
from chemflow2.units.mixer import Mixer
from chemflow2.units.reactor import Reactor
from chemflow2.units.separator import Separator
from chemflow2.units.splitter import Splitter

__all__ = ["Mixer", "Reactor", "Separator", "Splitter", "GibbsReactor", "CanteraError"]
