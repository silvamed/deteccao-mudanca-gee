// Extracao dos endmembers a partir da propria imagem
// PFC de Engenharia Cartografica, IME - Pedro Medeiros e Pablo Artte
//
// Cada endmember e a media dos pixels puros dentro de poligonos desenhados no
// mapa, sobre a mesma composicao do app (Sentinel-2 SR, Cloud Score+ >= 0,55,
// jun-set). Usamos a mediana das composicoes de 2019 a 2025, para ter um
// conjunto unico de endmembers para a serie inteira.
//
// Como usar:
//   1. Crie camadas de geometria com estes nomes:
//        vegetacao - mata densa e continua, em encosta iluminada
//        solo      - solo exposto dentro da propriedade (manchas de 30 x 30 m
//                    ou mais); nada de asfalto, cascalho ou telhado
//        sombra    - opcional (agua limpa e funda); sem ela fica 0,01
//      Varios poligonos pequenos por classe, por dentro das manchas.
//   2. Se quiser medir o erro de ajuste dentro do perimetro, preencha AOI_ASSET.
//   3. Rode e copie do console os valores para o PARAM do app.
//   Isso deve ser feito antes de sortear a amostra de validacao.

var ANO_INI = 2019, ANO_FIM = 2025;
var MES_INI = 6, MES_FIM = 9;
var BANDAS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];
var CS_LIMIAR = 0.55;
var ESCALA = 10;
var CRS = 'EPSG:31983';   // SIRGAS 2000 / UTM 23S (Rio de Janeiro)
var AOI_ASSET = '';       // opcional, ex.: 'projects/.../assets/limite_om'

// Endmembers da primeira versao do app, so para comparar o erro de ajuste.
var EM_ANTIGOS = {
  veg:    [0.030, 0.055, 0.030, 0.450, 0.160, 0.060],
  solo:   [0.110, 0.150, 0.220, 0.320, 0.420, 0.380],
  sombra: [0.010, 0.010, 0.010, 0.012, 0.010, 0.010]
};

var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var CS = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

if (typeof vegetacao === 'undefined' || typeof solo === 'undefined') {
  throw new Error('Desenhe as camadas "vegetacao" e "solo" antes de rodar.');
}
var temSombra = (typeof sombra !== 'undefined');

// a camada pode ser Geometry ou FeatureCollection
function geom(x) {
  return (x instanceof ee.FeatureCollection) ? x.geometry() : ee.Geometry(x);
}
var gVeg = geom(vegetacao);
var gSolo = geom(solo);
var gSombra = temSombra ? geom(sombra) : null;

var aoi = AOI_ASSET ? ee.FeatureCollection(AOI_ASSET).geometry()
                    : gVeg.union(gSolo, 1).bounds(1).buffer(2000, 100);
var regiao = aoi.buffer(2000, 100);

// Composicao de referencia (mesmo criterio do app)
function mascarar(image) {
  var img = ee.Image(image);
  return img.updateMask(img.select('cs_cdf').gte(CS_LIMIAR))
            .select(BANDAS).divide(10000);
}

var anuais = ee.ImageCollection(ee.List.sequence(ANO_INI, ANO_FIM).map(function(a) {
  var ini = ee.Date.fromYMD(a, MES_INI, 1);
  var fim = ee.Date.fromYMD(a, MES_FIM, 1).advance(1, 'month');
  return ee.ImageCollection(S2.filterDate(ini, fim).filterBounds(regiao)
                              .linkCollection(CS, ['cs_cdf'])
                              .map(mascarar)).median();
}));
var ref = anuais.median();
var ndvi = ref.normalizedDifference(['B8', 'B4']).rename('ndvi');

// Media, desvio padrao e numero de pixels de cada classe
function estatisticas(g) {
  return ref.addBands(ndvi).reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true)
                              .combine(ee.Reducer.count(), '', true),
    geometry: g, scale: ESCALA, crs: CRS, maxPixels: 1e9
  });
}

var pedido = {veg: estatisticas(gVeg), solo: estatisticas(gSolo)};
if (temSombra) { pedido.sombra = estatisticas(gSombra); }

// Erro medio quadratico do ajuste do MLME
function rmsMedio(em) {
  var f = ref.unmix([em.veg, em.solo, em.sombra], true, true);
  var soma = null;
  for (var i = 0; i < BANDAS.length; i++) {
    var modelo = f.select(0).multiply(em.veg[i])
      .add(f.select(1).multiply(em.solo[i]))
      .add(f.select(2).multiply(em.sombra[i]));
    var res = ref.select(BANDAS[i]).subtract(modelo);
    soma = (soma === null) ? res.multiply(res) : soma.add(res.multiply(res));
  }
  return soma.divide(BANDAS.length).sqrt().rename('rms').reduceRegion({
    reducer: ee.Reducer.mean(), geometry: aoi, scale: 30, crs: CRS,
    maxPixels: 1e10, bestEffort: true
  }).get('rms');
}

function tres(v) { return (Math.round(v * 1000) / 1000).toFixed(3); }
function linha(nome, lista) {
  return nome + '[' + lista.map(tres).join(', ') + '],';
}

ee.Dictionary(pedido).evaluate(function(r, err) {
  if (err) { print('Erro:', err); return; }

  function vetor(est) { return BANDAS.map(function(b) { return est[b + '_mean']; }); }
  function desvios(est) { return BANDAS.map(function(b) { return est[b + '_stdDev']; }); }

  var novos = {
    veg: vetor(r.veg),
    solo: vetor(r.solo),
    sombra: temSombra ? vetor(r.sombra) : EM_ANTIGOS.sombra
  };

  print('Pixels usados: vegetacao ' + r.veg.B4_count + ', solo ' + r.solo.B4_count +
        (temSombra ? ', sombra ' + r.sombra.B4_count : ', sombra fotometrica (0,01)'));
  print('NDVI medio: vegetacao ' + tres(r.veg.ndvi_mean) + ', solo ' + tres(r.solo.ndvi_mean));
  if (r.veg.ndvi_mean < 0.80) {
    print('ATENCAO: NDVI medio da vegetacao abaixo de 0,80. Os poligonos podem incluir ' +
          'sombra de encosta, borda ou vegetacao aberta; redesenhe.');
  }
  if (r.solo.ndvi_mean > 0.25) {
    print('ATENCAO: NDVI medio do solo acima de 0,25. Os poligonos incluem vegetacao; ' +
          'redesenhe dentro das manchas de solo exposto.');
  }
  print('Desvio padrao por banda (reflectancia):');
  print('  vegetacao ' + desvios(r.veg).map(tres).join(', '));
  print('  solo      ' + desvios(r.solo).map(tres).join(', '));

  print('--- Valores para colar no app (PARAM), ordem B2, B3, B4, B8, B11, B12 ---');
  print(linha('em_veg:    ', novos.veg));
  print(linha('em_solo:   ', novos.solo));
  print(linha('em_sombra: ', novos.sombra));

  ee.Dictionary({antigos: rmsMedio(EM_ANTIGOS), novos: rmsMedio(novos)})
    .evaluate(function(q, err2) {
      if (err2) { print('Erro no RMS:', err2); return; }
      print('RMS medio do ajuste: antigos ' + q.antigos.toFixed(4) +
            ' | novos ' + q.novos.toFixed(4) +
            (q.novos <= q.antigos ? '  (novos ajustam melhor ou igual)'
                                  : '  (novos ajustam PIOR: confira os poligonos)'));
    });

  // CSV com os valores obtidos
  var tabela = ee.FeatureCollection(['veg', 'solo', 'sombra'].map(function(k) {
    var props = {componente: k, origem: (k === 'sombra' && !temSombra) ? 'fotometrica' : 'imagem'};
    for (var i = 0; i < BANDAS.length; i++) { props[BANDAS[i]] = novos[k][i]; }
    return ee.Feature(null, props);
  }));
  Export.table.toDrive({collection: tabela, description: 'endmembers_imagem',
                        folder: 'PFC_degradacao', fileFormat: 'CSV'});
});

// Camadas para conferencia
Map.addLayer(ref, {bands: ['B4', 'B3', 'B2'], min: 0.02, max: 0.30}, 'Referencia 2019-2025 (jun-set)');
Map.addLayer(ndvi, {min: 0, max: 0.9, palette: ['b35806', 'f7f7f7', '1b7837']}, 'NDVI', false);
Map.addLayer(ee.FeatureCollection([ee.Feature(gVeg)]), {color: '00ff00'}, 'vegetacao');
Map.addLayer(ee.FeatureCollection([ee.Feature(gSolo)]), {color: 'ff8800'}, 'solo');
if (temSombra) { Map.addLayer(ee.FeatureCollection([ee.Feature(gSombra)]), {color: '0000ff'}, 'sombra'); }
gVeg.centroid(10).coordinates().evaluate(function(c) { if (c) { Map.setCenter(c[0], c[1], 13); } });
