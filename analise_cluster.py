import pandas as pd
import numpy as np
import re
import time
import os
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
# Adiciona os novos algoritmos para uma análise mais completa
from sklearn.cluster import (
    KMeans, 
    DBSCAN, 
    AgglomerativeClustering, 
    Birch, 
    SpectralClustering,
    MiniBatchKMeans  # <-- ADICIONADO O ALGORITMO RÁPIDO
)
# V-Measure é uma métrica externa (usa gabarito), similar ao F1-Score
from sklearn.metrics import (
    silhouette_score,      # Métrica Interna (não-supervisionada)
    davies_bouldin_score,  # Métrica Interna (não-supervisionada)
    v_measure_score        # Métrica Externa (supervisionada, nosso gabarito)
)
import warnings
import gc  # Garbage Collector para otimizar memória

# Suprimir avisos
warnings.filterwarnings('ignore')

# ----------------------------------------------------------------------
# 1. CONFIGURAÇÃO INICIAL
# ----------------------------------------------------------------------
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
K_MER_SIZE = 2
SKIP_X = 1  # K-mer 2x2 com skip 1
N_COMPONENTS = 300

# ----------------------------------------------------------------------
# 2. FUNÇÕES DE PROCESSAMENTO (Sem alterações)
# ----------------------------------------------------------------------
def parse_fasta(file_path):
    """Lê o arquivo FASTA e extrai a classe SCOPe (gabarito) e a sequência."""
    data = []
    current_header = None
    current_sequence = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_header and current_sequence:
                        match = re.search(r'>\w+\s+(\w\.\d+\.\d+\.\d+)', current_header)
                        if match:
                            scop_class = match.group(1)
                            main_class = scop_class[0]
                            data.append({
                                'scop_class_full': scop_class,
                                'main_class': main_class,
                                'sequence': "".join(current_sequence)
                            })
                    current_header = line
                    current_sequence = []
                elif line:
                    current_sequence.append(line.upper().replace(' ', ''))
        if current_header and current_sequence:
            match = re.search(r'>\w+\s+(\w\.\d+\.\d+\.\d+)', current_header)
            if match:
                scop_class = match.group(1)
                main_class = scop_class[0]
                data.append({
                    'scop_class_full': scop_class,
                    'main_class': main_class,
                    'sequence': "".join(current_sequence)
                })
        return pd.DataFrame(data)
    except FileNotFoundError:
        print(f"ERRO: Arquivo '{file_path}' não encontrado. Verifique o caminho.")
        return pd.DataFrame()
    except Exception as e:
        print(f"Ocorreu um erro durante o processamento do FASTA: {e}")
        return pd.DataFrame()

def generate_kmer_features(sequences, k=K_MER_SIZE, skip=SKIP_X):
    """Gera a matriz binária de presença/ausência de k-mers com skip."""
    kmer_list = [aa1 + aa2 for aa1 in AMINO_ACIDS for aa2 in AMINO_ACIDS]
    kmer_df = pd.DataFrame(0, index=range(len(sequences)), columns=kmer_list)
    print(f"  -> Total de {len(kmer_list)} atributos (k-mers) gerados.")
    
    for idx, seq in enumerate(sequences):
        clean_seq = ''.join(c for c in seq if c in AMINO_ACIDS)
        window_size = k + skip
        if len(clean_seq) < window_size:
            continue
        for i in range(len(clean_seq) - (window_size - 1)):
            kmer = clean_seq[i] + clean_seq[i + k - 1 + skip]
            if kmer in kmer_df.columns:
                kmer_df.loc[idx, kmer] = 1
    return kmer_df

def evaluate_clustering(estimator, X, y_true_encoded, algorithm_name):
    """Executa o algoritmo, mede o tempo e calcula as métricas."""
    start_time = time.time()
    try:
        y_pred = estimator.fit_predict(X)
    except Exception:
        return {
            'Algoritmo': algorithm_name, 'Nº Clusters': 0, 'Silhouette': -2,
            'Davies-Bouldin': 9999, 'V-Measure (Externo)': 0.0, 'Tempo (s)': time.time() - start_time
        }
    end_time = time.time()
    valid_indices = y_pred != -1
    X_valid, y_true_valid, y_pred_valid = X[valid_indices], y_true_encoded[valid_indices], y_pred[valid_indices]
    n_clusters = len(np.unique(y_pred_valid))
    if n_clusters < 2:
        sil_score, db_score, v_measure = -1, 999, 0.0
    else:
        sil_score = silhouette_score(X_valid, y_pred_valid)
        db_score = davies_bouldin_score(X_valid, y_pred_valid)
        v_measure = v_measure_score(y_true_valid, y_pred_valid)
    return {
        'Algoritmo': algorithm_name, 'Nº Clusters': n_clusters, 'Silhouette': sil_score,
        'Davies-Bouldin': db_score, 'V-Measure (Externo)': v_measure, 'Tempo (s)': end_time - start_time
    }

# ----------------------------------------------------------------------
# 3. EXECUÇÃO DO PROJETO
# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    # ETAPA 0: CAMINHOS DOS ARQUIVOS (Estão corretos, segundo seu terminal)
    FASTA_FILES = [
        r'G:\Outros computadores\JGFP NOTEBOOK\Documents\Faculdade\Topicos em Software\ExAnaliseDeDados\astral-scopedom-seqres-gd-sel-gs-bib-95-2.08.fa',
        r'G:\Outros computadores\JGFP NOTEBOOK\Documents\Faculdade\Topicos em Software\ExAnaliseDeDados\astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa'
    ]

    for fasta_path in FASTA_FILES:
        file_name = os.path.basename(fasta_path)
        print("\n" + "="*80)
        print(f"INICIANDO ANÁLISE PARA O ARQUIVO: {file_name}")
        print("="*80 + "\n")

        # ETAPA A: CARREGAMENTO
        print("--- ETAPA A: CARREGAMENTO DE DADOS ---")
        df_astral = parse_fasta(fasta_path)
        if df_astral.empty: continue
        print(f"-> Dados carregados: {len(df_astral)} sequências.")
        y_true = df_astral['main_class']
        num_classes = len(y_true.unique())
        print(f"-> Classes SCOPe (Gabarito): {num_classes} classes ({', '.join(y_true.unique())}).")
        if num_classes < 2:
            print(f"\nAVISO: Encontrada apenas {num_classes} classe. Pulando clusterização.")
            continue
        le = LabelEncoder()
        y_true_encoded = le.fit_transform(y_true)

        # ETAPA B: EXTRAÇÃO DE K-MER
        print("\n--- ETAPA B: EXTRAÇÃO DE ATRIBUTOS (K-MER BINÁRIO) ---")
        X_binario = generate_kmer_features(df_astral['sequence'])
        print(f"-> Matriz Binária gerada: {X_binario.shape}")

        # ETAPA C: REDUÇÃO DE DIMENSIONALIDADE (PCA)
        print("\n--- ETAPA C: REDUÇÃO DE DIMENSIONALIDADE (PCA) ---")
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_binario)
        pca = PCA(n_components=N_COMPONENTS)
        X = pca.fit_transform(X_scaled)
        explained_var = np.sum(pca.explained_variance_ratio_) * 100
        print(f"-> PCA finalizado. Matriz de atributos X: {X.shape}")
        print(f"-> Variância explicada com {N_COMPONENTS} componentes: {explained_var:.2f}%")
        del X_binario, X_scaled
        gc.collect()
        print("-> Matrizes intermediárias removidas da memória.")

        # --- ETAPA D MODIFICADA ---
        print("\n--- ETAPA D: TESTE DE ALGORITMOS (BASELINE) ---")
        
        # Lista de algoritmos que escalam para N=35k
        algorithms = [
            ('KMeans (K={})'.format(num_classes), KMeans(n_clusters=num_classes, random_state=42, n_init='auto')),
            ('MiniBatchKMeans (K={})'.format(num_classes), MiniBatchKMeans(n_clusters=num_classes, random_state=42, n_init='auto')),
            ('Birch (K={})'.format(num_classes), Birch(n_clusters=num_classes)),
            ('DBSCAN (Base)', DBSCAN(eps=5, min_samples=25, n_jobs=-1)), 
        ]
        # (Agglomerative e Spectral foram removidos por não escalarem)
        # --- FIM DA MODIFICAÇÃO ---

        all_results = []
        for name, estimator in algorithms:
            print(f"  -> Rodando: {name}...")
            results = evaluate_clustering(estimator, X, y_true_encoded, name)
            all_results.append(results)
        
        results_df = pd.DataFrame(all_results)
        print("\n--- RESULTADOS DA AVALIAÇÃO BASELINE ---")
        print(results_df.sort_values(by='V-Measure (Externo)', ascending=False).to_markdown(index=False))

        # ETAPA F: CORRELAÇÃO DE MÉTRICAS (A ESCOLHA NÃO-SUBJETIVA)
        print("\n" + "="*80)
        print("ETAPA F: CORRELAÇÃO DE MÉTRICAS (A ESCOLHA NÃO-SUBJETIVA)")
        print("="*80 + "\n")
        valid_results_df = results_df[results_df['Nº Clusters'] > 1].copy()
        if valid_results_df.empty:
            print("Nenhum algoritmo produziu clusters válidos. Impossível correlacionar.")
            continue
        metrics_to_correlate = valid_results_df[['Silhouette', 'Davies-Bouldin', 'V-Measure (Externo)']]
        correlation_matrix = metrics_to_correlate.corr()
        print("--- Matriz de Correlação (Internas vs Externa) ---")
        print(correlation_matrix.to_markdown())
        external_corr = correlation_matrix['V-Measure (Externo)'].drop('V-Measure (Externo)')
        print("\n--- Correlação com a Métrica Externa (V-Measure) ---")
        print(external_corr.to_markdown())
        best_internal_metric_name = external_corr.abs().idxmax()
        best_internal_metric_corr = external_corr.loc[best_internal_metric_name]
        print("\n--- CONCLUSÃO DA METODOLOGIA ---")
        print(f"A métrica interna (não-supervisionada) mais confiável é: {best_internal_metric_name}")
        print(f"Justificativa: Ela possui a maior correlação absoluta ({best_internal_metric_corr:.4f}) com a métrica externa (V-Measure).")
        print(f"A partir de agora, usaremos '{best_internal_metric_name}' como nossa 'bússola' para otimizar os parâmetros.")

        # ETAPA G: OTIMIZAÇÃO USANDO A MÉTRICA INTERNA VALIDADA
        print("\n--- ETAPA G: OTIMIZAÇÃO (ex: K-Means) USANDO A MÉTRICA VALIDADA ---")
        ascending_order = False
        if best_internal_metric_name == 'Davies-Bouldin':
            ascending_order = True
        K_RANGE = range(max(2, num_classes - 5), num_classes + 6)
        print(f"  -> Tunando K-Means com K variando de {K_RANGE.start} a {K_RANGE.stop-1}...")
        print(f"  -> Objetivo: Encontrar o 'k' com o melhor score de '{best_internal_metric_name}'")

        k_means_tuning_results = []
        for k in K_RANGE:
            if k < 2: continue
            km = KMeans(n_clusters=k, random_state=42, n_init='auto')
            results = evaluate_clustering(km, X, y_true_encoded, f'KMeans (k={k})')
            k_means_tuning_results.append(results)
        
        if k_means_tuning_results:
            tuning_df = pd.DataFrame(k_means_tuning_results)
            best_row = tuning_df.sort_values(by=best_internal_metric_name, ascending=ascending_order).iloc[0]
            print(f"\n--- MELHOR CONFIGURAÇÃO ENCONTRADA (Baseado em {best_internal_metric_name}) ---")
            print(best_row.to_markdown())
            best_k = best_row['Nº Clusters']
            best_internal_score = best_row[best_internal_metric_name]
            print(f"\nSugestão: O melhor K foi {best_k}, que atingiu o score de {best_internal_score:.4f} (usando {best_internal_metric_name}).")
        else:
            print("  -> Nenhum resultado de tuning do K-Means para exibir.")